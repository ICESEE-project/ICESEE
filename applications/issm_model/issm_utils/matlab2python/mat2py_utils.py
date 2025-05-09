# ====================================================================
# @author: Brian Kyanjo
# @description: Matlab-python server launcher for ISSM model and other helper functions
# @date: 2025-04-16
# ====================================================================

# --- Imports ---
import os
import sys
import subprocess
import numpy as np
import time
import signal
import psutil
import signal
import platform

class _MatlabServer:
    """A class to manage a MATLAB server for running ISSM models."""

    def __init__(self,color=0, matlab_path="matlab", cmdfile="cmdfile", statusfile="statusfile",verbose=False):
        """Initialize the MATLAB server configuration."""
        self.matlab_path = matlab_path
        self.cmdfile = os.path.abspath(f"{cmdfile}_{color}.txt")
        self.statusfile = os.path.abspath(f"{statusfile}_{color}.txt")
        self.process = None
        self.verbose = verbose

    # def kill_matlab_processes(self):
    #     matlab_count = 0
    #     for proc in psutil.process_iter(['name', 'pid']):
    #         try:
    #             # Check for MATLAB processes (name varies by OS)
    #             if 'matlab' in proc.info['name'].lower() or 'MATLAB' in proc.info['name']:
    #                 print(f"Found MATLAB process: {proc.info['name']} (PID: {proc.info['pid']})")
    #                 # Terminate the process
    #                 if platform.system() == "Windows":
    #                     proc.terminate()  # Windows uses terminate
    #                 else:
    #                     proc.send_signal(signal.SIGTERM)  # Unix uses SIGTERM
    #                 matlab_count += 1
    #         except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
    #             continue
    
    #     if matlab_count == 0:
    #         print("No MATLAB processes found running.")
    #     else:
    #         print(f"Terminated {matlab_count} MATLAB process(es).")


    def kill_matlab_processes(self):
        import psutil
        import platform
        import signal
        matlab_count = 0
        for proc in psutil.process_iter(['name', 'pid', 'cmdline']):
            try:
                # Check for MATLAB processes (name varies by OS)
                if 'matlab' in proc.info['name'].lower() or 'MATLAB' in proc.info['name']:
                    if self.verbose:
                        print(f"Found MATLAB process: {proc.info['name']} (PID: {proc.info['pid']})")
                    
                    # Get command-line arguments to check for GUI-related flags
                    cmdline = proc.info['cmdline']
                    is_gui = True  # Assume GUI unless proven otherwise
                    
                    # Check command-line arguments for non-GUI flags
                    if cmdline and any(flag in cmdline for flag in ['-nodisplay', '-nodesktop']):
                        is_gui = False  # Non-GUI instance
                    
                    if not is_gui:
                        # Terminate non-GUI process
                        if platform.system() == "Windows":
                            proc.terminate()  # Windows uses terminate
                        else:
                            proc.send_signal(signal.SIGTERM)  # Unix uses SIGTERM
                        matlab_count += 1

                        if self.verbose:
                            print(f"Terminated MATLAB process (PID: {proc.info['pid']})")
                    else:
                        print(f"Skipped GUI MATLAB process (PID: {proc.info['pid']})")
                        
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                continue
    
        if self.verbose:        
            if matlab_count == 0:
                print("No non-GUI MATLAB processes found to terminate.")
            else:
                print(f"Terminated {matlab_count} non-GUI MATLAB process(es).")


    def launch(self):
        """Launch MATLAB server and wait for it to be ready."""
        try:
            self.kill_matlab_processes()
        except Exception as e:
            print(f"An error occurred: {e}")
            
        if self.verbose:
            print("[Launcher] Starting MATLAB server...")
            print(f"[Launcher] Command file: {self.cmdfile}")
            print(f"[Launcher] Status file: {self.statusfile}")
            
        # Clean up old files if they exist
        for f in [self.cmdfile, self.statusfile]:
            if os.path.exists(f):
                os.remove(f)
        
        try:
            # Launch MATLAB in background with redirected I/O
            # matlab_cmd = f"{self.matlab_path} -nodesktop -nosplash -nojvm -r \"matlab_server('{self.cmdfile}', '{self.statusfile}, {int(self.verbose)}')\""
            matlab_cmd = f"{self.matlab_path} -nodesktop -nosplash -nojvm -r \"matlab_server('{self.cmdfile}', '{self.statusfile}')\""
            self.process = subprocess.Popen(
                matlab_cmd,
                shell=True,
                stdout=subprocess.PIPE,  # Redirect stdout
                stderr=subprocess.PIPE,  # Redirect stderr
                stdin=subprocess.PIPE,   # Redirect stdin
                preexec_fn=os.setsid    # Create new process group to handle signals
            )
            
            # Wait for server to signal readiness
            timeout = 10  # seconds
            start_time = time.time()
            while not os.path.exists(self.statusfile):
                if time.time() - start_time > timeout:
                    print("[Launcher] Error: MATLAB server failed to start within timeout.")
                    os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                    sys.exit(1)
                time.sleep(0.5)
            
            with open(self.statusfile, 'r') as f:
                status = f.read().strip()
            if status != 'ready':
                print(f"[Launcher] Error: Unexpected status '{status}'.")
                os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                sys.exit(1)
            
            if self.verbose:
                print("[Launcher] MATLAB server is ready.")
        except Exception as e:
            print(f"[Launcher] Error launching MATLAB server: {e}")
            if self.process:
                os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
            sys.exit(1)
            

    def send_command(self, command, timeout=300):
        """Send a command to MATLAB and wait for it to be processed."""
        if self.verbose:
            print(f"[Launcher] Sending command: {command}")
        with open(self.cmdfile, 'w') as f:
            f.write(command)
        
        # Wait for command to be processed (file deleted)
        start_time = time.time()
        sleep_time = 0.2 # start fast and slow down
        max_sleep  = 5.0 # don't sleep more than this
        while os.path.exists(self.cmdfile):
            elapsed_time = time.time() - start_time
            if elapsed_time > timeout:
                print("[Launcher] Error: Command execution timed out.")
                return False
            # if self.verbose:
            #     print("[Launcher] Waiting for command to be processed...")
            
            time.sleep(sleep_time)

            # adjust sleep time to slow down if needed up to max_sleep
            # sleep_time = min(sleep_time * 1.5, max_sleep)
            sleep_time = min(sleep_time + (elapsed_time / 10.0), max_sleep)
            if sleep_time == max_sleep:
                # step up the sleep time
                # sleep_time = 0.2
                print("[Launcher] Warning: Slow command processing detected.")
        
        print("[Launcher] Command processed successfully.")
        return True

    def shutdown(self):
        """Attempt to gracefully shut down the MATLAB server."""
        if self.verbose:
            print("[Launcher] Attempting to shut down MATLAB server...")

        if self.send_command("exit"):
            try:
                # Capture output to diagnose issues
                stdout, stderr = self.process.communicate(timeout=5)
                if self.verbose:
                    if stdout:
                        print("[MATLAB stdout]", stdout.decode())
                    if stderr:
                        print("[MATLAB stderr]", stderr.decode())
                    print("[Launcher] MATLAB server shut down successfully.")
            except subprocess.TimeoutExpired:
                print("[Launcher] Warning: MATLAB process did not terminate in time, forcing termination.")
                os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                self.process.wait(timeout=5)
                if self.verbose:
                    print("[Launcher] MATLAB server terminated.")
        else:
            print("[Launcher] Error: Failed to shut down MATLAB server gracefully.")
            os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
            self.process.wait(timeout=5)
            print("[Launcher] MATLAB server terminated.")
    
    def reset_terminal(self):
        """Reset terminal settings to restore normal behavior, if applicable."""
        if not sys.stdin.isatty() or any(key in os.environ for key in ["MPIEXEC", "OMPI_COMM_WORLD_RANK", "PMI_RANK", "SLURM_JOB_ID"]):
            if self.verbose:
                print("[Launcher] Skipping terminal reset (non-interactive or MPI environment).", file=sys.stderr, flush=True)
            return
        try:
            subprocess.run(['stty', 'sane'], check=True)
            if self.verbose:
                print("[Launcher] Terminal settings reset successfully.", file=sys.stderr, flush=True)
        except subprocess.CalledProcessError as e:
            print(f"[Launcher] Warning: Failed to reset terminal settings: {e}", file=sys.stderr, flush=True)

    # def reset_terminal(self):
    #     """Reset terminal settings to restore normal behavior."""
    #     try:
    #         subprocess.run(['stty', 'sane'], check=True)
    #         if self.verbose:
    #             print("[Launcher] Terminal settings reset successfully.")
    #     except subprocess.CalledProcessError:
    #         print("[Launcher] Warning: Failed to reset terminal settings.")

#  ---- end of MatlabServer class ----

# Lets use inheritance to create a new class that will manage the matlab server
class MatlabServer:
    """A class to manage a MATLAB server for running ISSM models both in parallel and serial."""
    def __init__(self, color=0, Nens=None, comm=None, matlab_path="matlab", cmdfile="cmdfile", statusfile="statusfile", verbose=False):
        """Initialize the MATLAB server configuration."""
        self.comm = comm
        self.rank = self.comm.Get_rank() if self.comm else 0
        # self.size = self.comm.Get_size() if self.comm else 1
        self.Nens = Nens
        self.color = color
        self.matlab_path = matlab_path
        self.cmdfile = cmdfile
        self.statusfile = statusfile
        self.verbose = verbose
        self._server = None

        from mpi4py import MPI
        self.size = MPI.COMM_WORLD.Get_size() if self.comm else 1

    def _check_conditions(self):
        """Check if conditions are met to access _MatlabServer."""
        if self.Nens is None:
            return False
        if self.Nens >= self.size:
            return True
        if self.Nens < self.size:
            raise ValueError(
                f"Nens ({self.Nens}) must be greater than or equal to the size of the MPI communicator ({self.size} set model_nprocs for the remaining resources if you want the coupled model to run in parallel). "
            )
        return False

    def __getattr__(self, name):
        """Delegate method calls to _MatlabServer if conditions are met."""
        if self._check_conditions():
            if self._server is None:
                self._server = _MatlabServer(
                    color=self.color,
                    matlab_path=self.matlab_path,
                    cmdfile=self.cmdfile,
                    statusfile=self.statusfile,
                    verbose=self.verbose
                )
                #  call the launch method to start the server
                self._server.launch()
            return getattr(self._server, name)
        # else:
        #     pass
        raise AttributeError(f"Method '{name}' not available: conditions not met (Nens={self.Nens}, size={self.size}, rank={self.rank})")

# --- Subprocess Command Runner ---

def subprocess_cmd_run(issm_cmd, nprocs: int, verbose: bool = True):
    try:
        process = subprocess.Popen(
            issm_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )

        stdout, stderr = process.communicate()

        if verbose:
            stdout_lines = stdout.splitlines()
            trimmed_stdout = "\n".join(stdout_lines[9:])
            print(f"\n[ICESEE] ➤ Running ISSM with {nprocs} processors")
            print("------ ICESEE<->MATLAB STDOUT ------")
            print(trimmed_stdout.strip())

            if stderr.strip():
                print("------ ICESEE<->MATLAB STDERR ------")
                print(stderr.strip())

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, issm_cmd)

    except FileNotFoundError:
        print("❌ Error: MATLAB not found in PATH.")
    except subprocess.CalledProcessError as e:
        print(f"❌ MATLAB exited with error code {e.returncode}")
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        
# def subprocess_cmd_run(issm_cmd, nprocs: int, verbose: bool = True):
#     """
#     Run ISSM using a MATLAB script via subprocess.Popen.

#     Parameters:
#     - issm_cmd: Full command string to run ISSM in MATLAB
#     - nprocs: Number of processors to pass to runme.m (for display/debug)
#     - verbose: If True, print stdout (trimmed) and stderr (only if non-empty)
#     """
#     try:
#         process = subprocess.Popen(
#             issm_cmd,
#             shell=True,
#             stdout=subprocess.PIPE,
#             stderr=subprocess.PIPE,
#             universal_newlines=True
#         )

#         stdout, stderr = process.communicate()

#         if verbose:
#             stdout_lines = stdout.splitlines()
#             trimmed_stdout = "\n".join(stdout_lines[9:])  # Skip banner
#             print(f"\n[ICESEE] ➤ Running ISSM with {nprocs} processors")
#             print("------ ICESEE<->MATLAB STDOUT ------")
#             print(trimmed_stdout.strip())

#             if stderr.strip():  # Only print stderr if there's content
#                 print("------ ICESEE<->MATLAB STDERR ------")
#                 print(stderr.strip())

#         if process.returncode != 0:
#             # raise subprocess.CalledProcessError(
#             #     process.returncode, issm_cmd, output=stdout, stderr=stderr
#             # )
#             raise subprocess.CalledProcessError(process.returncode, issm_cmd)

#     except FileNotFoundError:
#         print("❌ Error: MATLAB not found in PATH.")
#     except subprocess.CalledProcessError as e:
#         print(f"❌ MATLAB exited with error code {e.returncode}")
#         if e.stderr.strip():
#             print("------ MATLAB STDERR ------")
#             print(e.stderr.strip())
#     except Exception as e:
#         print(f"❌ Unexpected error: {e}")
        
#  --- Add ISSM_DIR to sys.path ---
def add_issm_dir_to_sys_path(issm_dir=None):
    """
    Add ISSM_DIR and its subdirectories to sys.path.

    Parameters:
    - issm_dir: str or None
        The ISSM directory path. If None, it tries to get from environment variable 'ISSM_DIR'.
    """

    import os
    import sys

    if issm_dir is None:
        issm_dir = os.environ.get('ISSM_DIR')

    if not issm_dir:
        raise EnvironmentError("ISSM_DIR is not set. Please set the ISSM_DIR environment variable.")

    if not os.path.isdir(issm_dir):
        raise FileNotFoundError(f"The ISSM_DIR directory does not exist: {issm_dir}")

    for root, dirs, _ in os.walk(issm_dir):
        sys.path.insert(0, root)

    print(f"[ICESEE] Added ISSM directory and subdirectories from path: {issm_dir}")


# --- MATLAB Engine Initialization ---
# MATLAB Engine Initialization
def initialize_matlab_engine():
    """
    Initializes the MATLAB Engine for Python.

    Returns:
        eng: The MATLAB Engine instance.
    """
    try:
        import matlab.engine
        print("Starting MATLAB engine...")
        # Start a headless MATLAB engine without GUI
        eng = matlab.engine.start_matlab("-nodisplay -nosplash -nodesktop -nojvm")
        print("MATLAB engine started successfully.")
        return eng
    except ImportError:
        print("MATLAB Engine API for Python not found. Attempting to install...")

        try:
            # Find MATLAB root
            matlab_root = find_matlab_root()

            # Install the MATLAB Engine API for Python
            install_matlab_engine(matlab_root)

            # Retry importing and starting the MATLAB Engine
            import matlab.engine
            eng = matlab.engine.start_matlab("-nodisplay -nosplash -nodesktop -nojvm")
            print("MATLAB engine started successfully after installation.")
            return eng
        except Exception as e:
            raise RuntimeError(
                f"Failed to initialize MATLAB Engine API for Python: {e}\n"
                "Ensure MATLAB is installed, and its bin directory is added to your PATH.\n"
                "For instructions on installing the MATLAB Engine API for Python, see "
                "the official MATLAB documentation or the provided README.md."
            )

# # --- MATLAB Root Finder ---
def find_matlab_root():
    """
    Finds the MATLAB root directory by invoking MATLAB from the terminal.

    Returns:
        matlab_root (str): The root directory of the MATLAB installation.
    """
    try:
        # Run MATLAB in terminal mode to get matlabroot
        result = subprocess.run(
            ["matlab", "-batch", "disp(matlabroot)"],  # Run MATLAB with -batch mode
            text=True,
            capture_output=True,
            check=True
        )
        
        # Extract and clean the output
        matlab_root = result.stdout.strip()
        print(f"MATLAB root directory: {matlab_root}")
        return matlab_root
    except FileNotFoundError:
        print(
            "MATLAB is not available in the system's PATH. "
            "Ensure MATLAB is installed and its bin directory is in the PATH."
        )
    except subprocess.CalledProcessError as e:
        print(f"Error while executing MATLAB: {e.stderr.strip()}")
        raise

# --- MATLAB Engine Installation ---
def install_matlab_engine(matlab_root):
    """
    Installs the MATLAB Engine API for Python using the MATLAB root directory.

    Args:
        matlab_root (str): The root directory of the MATLAB installation.
    """
    try:
        # Save the current working directory
        current_dir = os.getcwd()

        # Path to the setup.py script for MATLAB Engine API for Python
        setup_path = os.path.join(matlab_root, "extern", "engines", "python")
        assert os.path.exists(setup_path), f"Setup path does not exist: {setup_path}"

        # Change to the setup directory
        os.chdir(setup_path)

        # Run the setup.py script to install the MATLAB Engine API
        print("Installing MATLAB Engine API for Python...")
        result = subprocess.run(
            ["python", "setup.py", "install", "--user"],
            text=True,
            capture_output=True,
            check=True
        )

        # Export the build directory to PYTHONPATH
        home_path = os.path.expanduser("~/")  # Adjust if needed
        os.environ["PYTHONPATH"] = f"{home_path}/lib:{os.environ.get('PYTHONPATH', '')}"

        print("MATLAB Engine API for Python installed successfully.")
    except AssertionError as e:
        print(f"AssertionError: {e}")
        raise
    except subprocess.CalledProcessError as e:
        print(f"Error while installing MATLAB Engine API for Python:\n{e.stderr.strip()}")
        raise
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise
    finally:
        # Change back to the original directory
        os.chdir(current_dir)

# Example usage
if __name__ == "__main__":
    matlab_root = find_matlab_root()
    # install_matlab_engine(matlab_root)
