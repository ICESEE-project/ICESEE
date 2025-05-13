# ====================================================================
# @author: Brian Kyanjo
# @description: Matlab-python server launcher utilities for ISSM model and helper functions
# @date: 2025-04-16
# ====================================================================

import traceback
import sys
import os
import logging
import signal

# Global server reference and shutdown state
_global_server = None
_server_shutdown = False
_signal_handler_set = False  # Track if signal handler is set

def setup_server_shutdown(server, comm, verbose=False):
    """Set up global server reference and SIGINT handler for cleanup on KeyboardInterrupt.

    Configures a signal handler for SIGINT to gracefully shut down the MATLAB server.
    Only rank zero prints messages to stderr when verbose is True, but all ranks log to a file.

    Args:
        server: The MATLAB server instance to manage.
        comm (mpi4py.MPI.Comm or None): MPI communicator for rank checking. If None, assumes single process.
        verbose (bool, optional): If True, enables printing to stderr for rank zero. Defaults to False.
    """
    global _global_server, _server_shutdown, _signal_handler_set
    if _signal_handler_set:
        logging.debug("Global SIGINT handler already set up")
        return

    _global_server = server
    _server_shutdown = False

    # Configure logging to a file
    if os.path.exists("launcher_log.txt"):
        os.remove("launcher_log.txt")
    logging.basicConfig(
        filename="launcher_log.txt",
        level=logging.DEBUG,
        format="%(asctime)s [Launcher] %(message)s",
        force=True
    )

    def log_message(message):
        """Log a message to file and print to stderr if rank is 0 and verbose is True."""
        rank = comm.Get_rank() if comm is not None else 0
        logging.debug(message)  # Always log to file
        if verbose and rank == 0:
            print(f"[Launcher] {message}", file=sys.stderr, flush=True)

    # Signal handler for SIGINT
    def signal_handler(sig, frame):
        global _global_server, _server_shutdown
        log_message("Received SIGINT (KeyboardInterrupt) globally")
        if _global_server is not None and not _server_shutdown:
            log_message("Attempting to shut down MATLAB server globally")
            try:
                _global_server.shutdown()
                _server_shutdown = True
                # Skip reset_terminal in non-interactive or MPI environments
                mpi_env_vars = ["MPIEXEC", "OMPI_COMM_WORLD_RANK", "PMI_RANK", "SLURM_JOB_ID"]
                if sys.stdin.isatty() and not any(key in os.environ for key in mpi_env_vars):
                    log_message("Calling reset_terminal globally")
                    _global_server.reset_terminal()
                else:
                    log_message("Skipping reset_terminal globally (non-interactive or MPI environment)")
            except Exception as shutdown_error:
                log_message(f"Failed to shutdown MATLAB server globally: {shutdown_error}")
            _global_server = None  # Clear reference after shutdown
        log_message("Exiting due to SIGINT")
        sys.exit(1)

    signal.signal(signal.SIGINT, signal_handler)
    _signal_handler_set = True
    log_message("Global SIGINT handler set up for server shutdown")

def run_icesee_with_server(callable_func, server, shut_down_server=False, comm=None, verbose=True):
    """Run the specified callable with the given server, handling errors and cleanup.

    Executes a callable function, manages errors, and optionally shuts down the server.
    Only rank zero prints messages to stderr when verbose is True, but all ranks log to a file.

    Args:
        callable_func: A callable function or lambda to execute with the server.
        server: The MATLAB server instance to interact with.
        shut_down_server (bool, optional): If True, shuts down the server after execution. Defaults to False.
        comm (mpi4py.MPI.Comm or None, optional): MPI communicator for rank checking. Defaults to None.
        verbose (bool, optional): If True, enables printing to stderr for rank zero. Defaults to True.

    Returns:
        The result of the callable if successful, False if an error occurs.

    Usage:
        - For a MATLAB server method: run_icesee_with_server(lambda: server.send_command('your_command'), server)
        - For a regular function: run_icesee_with_server(function(args), server)
    """
    global _server_shutdown
    # Ensure logging is configured
    logging.basicConfig(
        filename="launcher_log.txt",
        level=logging.DEBUG,
        format="%(asctime)s [Launcher] %(message)s",
        force=True
    )

    def log_message(message):
        """Log a message to file and print to stderr if rank is 0 and verbose is True."""
        rank = comm.Get_rank() if comm is not None else 0
        logging.debug(message)  # Always log to file
        if verbose and rank == 0:
            print(f"[Launcher] {message}", file=sys.stderr, flush=True)

    # Ensure Python output is unbuffered
    os.environ["PYTHONUNBUFFERED"] = "1"
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', buffering=1)

    # Get function name for logging
    func_name = getattr(callable_func, '__name__', 'anonymous_function')
    
    success = False
    try:
        log_message(f"Starting {func_name}")
        result = callable_func()
        log_message(f"{func_name} result: {result}")
        if not result:
            raise RuntimeError(f"{func_name} failed.")
        success = True
        return result  # Return the callable's result

    except RuntimeError as e:
        log_message(f"RuntimeError: {e}")
        log_message(f"Full traceback: {traceback.format_exc()}")
        return False
    except Exception as e:
        log_message(f"Unexpected error: {e}")
        log_message(f"Full traceback: {traceback.format_exc()}")
        return False
    except KeyboardInterrupt:
        log_message("Process interrupted by user within run_icesee_with_server")
        return False
    finally:
        # Determine if shutdown is required
        should_shutdown = not (success and not shut_down_server)
        if verbose:
            log_message(f"Shutdown decision: should_shutdown={should_shutdown}, success={success}, shut_down_server={shut_down_server}")

        if should_shutdown:
            log_message("Attempting to shut down MATLAB server within run_icesee_with_server")
            try:
                server.shutdown()
                _server_shutdown = True
                # Skip reset_terminal in non-interactive or MPI environments
                mpi_env_vars = ["MPIEXEC", "OMPI_COMM_WORLD_RANK", "PMI_RANK", "SLURM_JOB_ID"]
                if sys.stdin.isatty() and not any(key in os.environ for key in mpi_env_vars):
                    log_message("Calling reset_terminal within run_icesee_with_server")
                    server.reset_terminal()
                else:
                    log_message("Skipping reset_terminal within run_icesee_with_server (non-interactive or MPI environment)")
            except Exception as shutdown_error:
                log_message(f"Failed to shutdown MATLAB server within run_icesee_with_server: {shutdown_error}")
        log_message("Exiting run_icesee_with_server")