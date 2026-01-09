#include <mpi.h>                                                                
#include <stdio.h>                                                              
#include <unistd.h>                                                             
#include <limits.h>                                                             
#include <stdlib.h>                                                             
#include <sys/syscall.h>                                                        
#include <sys/types.h>                                                          

int main(int argc, char** argv) {                                               
  // Initialize the MPI environment. The two arguments to MPI Init are not      
  // currently used by MPI implementations, but are there in case future        
  // implementations might need the arguments.                                  
  MPI_Init(NULL, NULL);                                                         

  // Get the number of processes                                                
  int world_size;                                                               
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);                                   

  // Get the rank of the process                                                
  int world_rank;                                                               
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);                                   

  // Get the name of the processor                                              
  char processor_name[MPI_MAX_PROCESSOR_NAME];                                  
  int name_len;                                                                 
  MPI_Get_processor_name(processor_name, &name_len);                            

  // Get CPU and NUMA information                                               
  unsigned int cpu, node;                                                       
  syscall(SYS_getcpu, &cpu, &node, NULL);                                       

  // Get mount namespace for container identity                                 
  char mnt_ns_buff[PATH_MAX];                                                   
  ssize_t mnt_ns_len = readlink("/proc/self/ns/mnt", mnt_ns_buff, sizeof(mnt_ns_buff)-1);
  if (mnt_ns_len == -1) {                                                       
    fprintf(stderr, "error getting mount namespace\n");                         
    exit(-1);                                                                   
  }                                                                             
  mnt_ns_buff[mnt_ns_len] = '\0';                                               

  // Print off a hello world message                                            
  printf("Hello world! Processor %s, Rank %d of %d, CPU %d, NUMA node %d, Namespace %s\n",
         processor_name, world_rank, world_size, cpu, node, mnt_ns_buff);       

  // Finalize the MPI environment. No more MPI calls can be made after this     
  MPI_Finalize();                                                               
}