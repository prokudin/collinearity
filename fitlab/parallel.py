from collections import namedtuple, deque
from multiprocessing import Process, current_process
from timeit import default_timer as timer
import logging
import psutil
import zmq
import cPickle as pickle
import numpy as np
import socket
import time
import os
from tools.config import conf
try:
  import mkl
except:
  pass

use_pickles = False

ZMQ_PULL = zmq.PULL  # @UndefinedVariable
ZMQ_PUSH = zmq.PUSH  # @UndefinedVariable
ZMQ_REQ = zmq.REQ  # @UndefinedVariable
ZMQ_REP = zmq.REP  # @UndefinedVariable
ZMQ_ROUTER = zmq.ROUTER  # @UndefinedVariable
ZMQ_POLLIN = zmq.POLLIN  # @UndefinedVariable
ZMQ_ROUTER_MANDATORY = zmq.ROUTER_MANDATORY  # @UndefinedVariable

def Server(ip):
    return ZmqServer(ip)

def Worker(ip):
    return ZmqWorker(ip)

def Broker(partition=None):
    return ZmqBroker(partition)

class ZmqBroker(object):

    def __init__(self, partition):
        if not partition:
            partition = psutil.cpu_count(logical=True)
        self.partition = partition
        self.process = None

    def run(self):
        context = zmq.Context(2)

        # Socket facing clients
        frontend = context.socket(ZMQ_ROUTER)
        frontend.bind("tcp://*:55550")

        # Socket facing services
        backend = context.socket(ZMQ_ROUTER)
        backend.bind("tcp://*:55551")
        backend.router_mandatory = 1
        #backend.setsockopt(ZMQ_ROUTER_MANDATORY, 1)

        # Idle workers
        ready_workers = deque()

        # FIFO Queue of unassigned work
        pending_work = deque()

        # Accumulation of results from workers, per-client
        results = {}

        poller = zmq.Poller()
        poller.register(backend, ZMQ_POLLIN)
        poller.register(frontend, ZMQ_POLLIN)

        # Used for debugging
        serial = 0
        workers = {}

        while True:
            sockets = dict(poller.poll())

            # Incoming from worker
            if backend in sockets:
                request = backend.recv_multipart()
                worker, _, client = request[:3]
                ready_workers.append(worker)

                # If this is a result message as opposed to a ready message, it will
                # contain more than three parts
                if len(request) > 3:
                    _, sn, result = request[3:]
                    results[client].append(result)

                    # If all workers have generated results, we can pass them to the client
                    if len(results[client]) == self.partition:
                        msg = [client, '']
                        msg.extend(results[client])
                        frontend.send_multipart(msg)

                else:
                    workers[worker] = client

            # Incoming from master
            if frontend in sockets:
                client, _, request = frontend.recv_multipart()
                results[client] = []
                pending_work.append([client, request, 0])

            # Can we assign some work?
            while pending_work and ready_workers:
                worker = ready_workers.popleft()
                client, request, offset = pending_work[0]
                serial += 1
                try:
                    backend.send_multipart([worker, '', client, '', repr(serial), repr(
                        offset), repr(self.partition), request], copy=True)
                except zmq.ZMQError as e:
                    continue
                offset += 1
                if offset < self.partition:
                    pending_work[0][2] = offset
                else:
                    pending_work.popleft()

    def run_subprocess(self):
        def run_runner():
            self.run()
        self.process = Process(target=run_runner, name='broker')
        self.process.start()

    def stop(self):
        if self.process:
            self.process.terminate()
            self.process.join(3)

class ZmqServer(object):

    def __init__(self, ip):
        self.broker_ip = ip
        self.theory = np.ndarray(0)
        self.received = False
        self.context = None
        self.worksock = None
        self.mproc_times = []
        self.param_times = []
        self.load_times = []
        self.store_times = []
        self.recv_times = []

    def assign_work(self):
        if not self.context:
            self.context = zmq.Context()
            self.worksock = self.context.socket(ZMQ_REQ)
            self.worksock.connect('tcp://' + self.broker_ip + ':55550')

        t = timer()
        state_tuples = []
        for x in ['pdf','gk','transversity','ffpi','ffk', 'collinspi','collinsk']:
            if x in conf:
                state_tuples.append((x, conf[x].get_state()))
        state = pickle.dumps(state_tuples, pickle.HIGHEST_PROTOCOL)
        self.worksock.send(state)
        self.received = False
        self.param_times.append(timer() - t)

    def wrap_mproc(self, obs_name, mproc):
        # Remember the result table offset of this observable
        first = self.theory.size
        self.theory = np.resize(self.theory, first + len(mproc.data))

        def run_mproc():
            a = timer()
            if not self.received:
                t = timer()
                data = self.worksock.recv_multipart()
                self.recv_times.append(timer() - t)
                t = timer()
                for chunk in data:
                    # Each chunk is an array of floats:
                    # chunk[0] is offset
                    # chunk[1] is stride
                    # chunk[2+n] is theory[offset+stride*n]
                    a = np.frombuffer(chunk, dtype=np.float64)
                    offset, stride = int(a[0]), int(a[1])
                    # print offset,stride,a.size
                    self.theory[offset::stride] = a[2:]
                self.received = True
                self.load_times.append(timer() - t)
            t = timer()
            tuples = [(mproc.data[i][0], mproc.data[i][1], self.theory[first + i])
                      for i in range(len(mproc.data))]
            self.store_times.append(timer() - t)
            self.mproc_times.append(timer() - a)

            #verify = mproc.singlecore()
            #verify.sort(key = lambda x: x[1])
            #tuples.sort(key = lambda x: x[1])
            # for i in range(len(verify)):
            #    if verify[i] != tuples[i]:
            #        logging.warn("Mismatch parallel [%d] %s != %s",i,verify[i],tuples[i])
            return tuples

        Wrapper = namedtuple('Wrapper', 'run')
        return Wrapper(run_mproc)

    def finis(self):
        pass

class ZmqWorker(object):

    def __init__(self, ip):
        self.ip = ip
        self.bigtable = []
        self.process = None
        self.workers = []

    def add_mproc(self, name, mproc):
        Row = namedtuple('Row', 'func entry')
        for entry in mproc.data:
            self.bigtable.append(Row(mproc.func, entry))

    def toil(self):
        mkl.set_num_threads(1)
        ctx = zmq.Context()
        worksock = ctx.socket(ZMQ_REQ)
        worksock.connect('tcp://' + self.ip + ':55551')
        worksock.send(socket.gethostname().split('.')
                      [0] + ':' + current_process().name)

        while True:
            address, _, serial, offset, stride, blob = worksock.recv_multipart()
            serial = int(serial)
            if stride == None:
                break
            offset, stride = int(offset), int(stride)
            state_tuples = pickle.loads(blob)
            for item in state_tuples:
                name, tuple = item[:]
                # logging.debug('name=%s',name)
                conf[name].set_state(tuple)

            results = [offset, stride]
            results.extend([row.func(row.entry)[2]
                            for row in self.bigtable[offset::stride]])
            r = np.array(results, dtype=np.float64)
            worksock.send_multipart([address, '', repr(serial), r])


    def run(self):

        if 'nprocs' in conf:
            nprocs = conf['nprocs']
        else:
            nprocs = psutil.cpu_count(logical=True)

        for i in range(nprocs):
            worker = Process(target=self.toil, name='w' + str(i))
            worker.start()
            psutil.Process(worker.ident).cpu_affinity([i])
            self.workers.append(worker)

        while len(self.workers):
            for w in list(self.workers):
                if not w.is_alive():
                    self.workers.remove(w)
            time.sleep(5)

    def run_subprocess(self):
        def run_runner():
            self.run()
        self.process = Process(target=run_runner, name='worker')
        self.process.start()

    def stop(self):
        for w in self.workers:
            w.terminate()
        for w in self.workers:
            w.join(3)
        if self.process:
            self.process.terminate()
            self.process.join(3)
