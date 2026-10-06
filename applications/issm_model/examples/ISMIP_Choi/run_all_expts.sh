mpirun -np 8 python run_da_issm.py --Nens=40 --model_nprocs=1 -F rebutal_experiments/param_ibf.yaml && \
sleep 5m && \
mpirun -np 8 python run_da_issm.py --Nens=40 --model_nprocs=1 -F rebutal_experiments/param_wbf.yaml && \
sleep 5m && \
mpirun -np 8 python run_da_issm.py --Nens=40 --model_nprocs=1 -F rebutal_experiments/param_ebf.yaml