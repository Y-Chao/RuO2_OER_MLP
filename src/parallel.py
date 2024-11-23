#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'

"""
This module is used to implement parallel computing.
Examples:
---------
If you want to run a light vasp task on 16/32 cores, you can use job scheduler to apply a 128 cores.
Then this module can help you to run 8/4 tasks on each core parallel. 
"""
import os
import traceback
from collections import Counter
import subprocess
import multiprocessing as mp

class Parallel:
    def __init__(self, ncores_per_task, hpc='NUS', workdir='.'):
        """
        Parameters:
        -----------
        ncores_per_task: int
            The number of cores per task.
        hpc: str
            The type of hpc, default is NUS.
        workdir: str
            The working directory, default is current directory.
        """
        self.workdir = workdir
        self.hpc = hpc
        self.ncores_per_task = ncores_per_task
    
    def set_process(self):
        """
        Set the allocated cores number.
        """
        # get the allocated nodes and cores by the job scheduler
        node_cores = self.get_process(self.hpc)
        # get the re-allocated chunk list
        chunk_list, total_processors = self.chunks(node_cores, self.ncores_per_task)
        print(" The chunk list of total processors: ", chunk_list)
        pool_size = len(chunk_list)
        print(" The processor pools number: ", pool_size)
        self.hostInfo = self.dump_hostfile(chunk_list)
        self.dump_hostfile(chunk_list)
        return self.hostInfo, pool_size, total_processors
    
    def get_process(self):
        """
        Get the allocated cores number.
        """
        node_cores = {}
        if self.hpc in ['NUS', 'NSCC']:
            hostfile = os.environ['PBS_NODEFILE']
            with open(hostfile, 'r') as f:
                nodes = f.read().splitlines()
                node_counts = Counter(nodes)
                node_cores = dict(node_counts)
            return node_cores
        else:
            raise ValueError('The hpc type is not supported.')
        
    def chunks(self, node_cores, n):
        """
        Yield successive n-sized chunks from node_cores.
        """
        allocated_processor = [ [node]*ncore_per_node for node, ncore_per_node in node_cores.items() ]
        total_processors = len(allocated_processor)
        print(" Available processors number: ", total_processors)
        # chunk the allocated processors
        chunk_list = [total_processors[i: i+n] for i in range(0, len(total_processors), n)]
        return chunk_list, total_processors
    
    def dump_hostfile(self, chunk_list):
        """
        Dump the hostfile according to the chunk list.
        Parameters:
        -----------
        chunk_list: list
            The list of chunk accroding to the cpu_per_jobs.
        Returns:
        --------
        hostInfo: list
            The list of (hostfile, core numebrs).
        """
        os.chdir(self.workdir)
        hostInfo = []
        # dump the hostfile accroding to the chunk list
        for i, chunk in enumerate(chunk_list):
            with open('.hostfile_%03d'%i, 'w') as f:
                for node in chunk:
                    f.write(node + '\n')
        hostInfo.append(('.hostfile_%03d'%i, len(chunk)))
        os.system('cat .hostfile_* > .hostfile')
        return hostInfo


"""
Then, we define some run functions.
"""
def run_prog_local(workdir, prog, ncores_per_task, hpc='NUS'):
    """
    Run the program on local machine.
    """
    try:
        os.chdir(workdir)
        mpi_prog = "mpirun -np %d " %ncores_per_task + prog
        with open('output', 'w') as f:
            subprocess.call(mpi_prog, stdout=f, stderr=f, shell=True, executable='/bin/bash')
        return
    except Exception as e:
        traceback.print_exc()
        raise e
    
def run_prog_hpc(work_dir, prog, hostInfo, root_dir='./', poolcount=0):
    """
    Run the program on hpc.
    The run_prog_hpc function should be called in the pool.apply_async.
    Examples:
    ---------
        task = 100
        exit = False
        pool = mp.Pool(pool_size)
        for i in range(len(task)):
            pool.apply_async(run_prog_hpc, args=(work_dir=i, prog=10, hostInfo=hostInfo, root_dir='./', poolcount=i))
        while not exit:
            time.sleep(10)
            if len(pool._cache) == 0:
                exit = True
                pool.terminate() # terminate all the processes in the pool, no active task will be processed
        pool.close() # close the pool, no more tasks can be added to the pool
        pool.join() # wait for all the tasks to be completed, usually called after close()
    """
    try:
        os.chdir(root_dir)
        cwd = os.getcwd()
        run_chuck = hostInfo[int(mp.current_process().name.split('-')[-1]) -1 - poolcount]
        os.chdir(work_dir)
        mf = os.path.join(cwd, run_chuck[0])

        mpi_prog = "mpirun -r ssh -machinefile %s -np %d" %(mf, run_chuck[1]) + prog
        with open('proginfo', 'w') as f1:
            f1.write('Current process: ' + mp.current_process().name + '\n')
            f1.write('Host file: ' + mf + '\n')
            f1.write('Run command: ' + mpi_prog + '\n')
        with open ('output', 'w') as fout:
            subprocess.call(mpi_prog, stdout=fout, stderr=fout, shell=True, executable='/bin/bash')
        return
    except Exception as e:
        traceback.print_exc()
        raise e

                
        



        
        
    
        
    
    
        



                