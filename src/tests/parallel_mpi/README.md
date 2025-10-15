# Parallel MPI Tests

Test suite for MPI parallelization functionality in ICESEE.

## Purpose

Validate parallel execution, communication patterns, and scalability of the MPI implementation.

## Test Coverage

Tests include:
- MPI initialization and finalization
- Ensemble distribution across processors
- Parallel I/O operations
- Communication patterns (scatter, gather, broadcast)
- Parallel forecast step
- Parallel analysis step
- Load balancing
- Scaling performance

## Running Tests

### Serial Test (for debugging)
```bash
python test_mpi_basic.py
```

### Parallel Tests
```bash
# 4 processes
mpirun -np 4 python -m pytest test_parallel_*.py

# 8 processes
mpirun -np 8 python -m pytest test_parallel_*.py
```

## Test Structure

MPI tests have special requirements:
- Must be run with mpirun/mpiexec
- Each process executes the test
- Assertions must be coordinated
- Cleanup is critical

## Performance Tests

### Strong Scaling
```bash
# Run with increasing processor counts
for np in 2 4 8 16; do
    mpirun -np $np python test_scaling.py --mode strong
done
```

### Weak Scaling
```bash
# Run with problem size scaling with processors
for np in 2 4 8 16; do
    mpirun -np $np python test_scaling.py --mode weak --size $((np*100))
done
```

## Common Issues

### Deadlocks
If tests hang, check for:
- Unmatched send/receive pairs
- Missing barriers
- Incorrect collective operations

### Load Imbalance
Monitor processor utilization:
```bash
mpirun -np 8 python test_load_balance.py --profile
```

## Requirements

- MPI implementation (OpenMPI, MPICH, Intel MPI)
- mpi4py
- pytest
- ICESEE core modules

## Expected Results

- All processors should complete successfully
- Communication overhead should be reasonable
- Near-linear speedup for embarrassingly parallel operations
- Acceptable efficiency for analysis step

## Debugging MPI Tests

### Verbose Output
```bash
mpirun -np 4 python -m pytest -v -s
```

### Single Process Debug
```bash
# Debug rank 0 only
mpirun -np 4 xterm -e gdb python test_parallel.py
```

### MPI Profiling
```bash
mpirun -np 8 --profile python test_performance.py
```

## Adding Tests

When adding MPI tests:
1. Ensure all processes participate
2. Use collective operations correctly
3. Validate results on all ranks
4. Clean up MPI resources
5. Handle edge cases (1 processor, many processors)

## Best Practices

- Test with different processor counts
- Verify correctness before performance
- Use timeouts to catch deadlocks
- Profile to identify bottlenecks
- Test both small and large problems
