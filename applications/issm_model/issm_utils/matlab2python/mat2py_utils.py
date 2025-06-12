# ====================================================================
# @author: Brian Kyanjo
# @description: Matlab-python server launcher for ISSM model and other helper functions
# @date: 2025-04-16
# ====================================================================


import os
import time
import subprocess
import sys
import signal
import psutil
import platform
import threading
import queue

class _MatlabServer:
    """A class to manage a MATLAB server for running ISSM models.

    This class provides functionality to launch a MATLAB server in non-GUI mode, send commands
    to it, and capture its output in real-time when verbose mode is enabled. In MPI environments,
    only the process with rank zero will print to the screen.

    Attributes:
        matlab_path (str): Path to the MATLAB executable.
        cmdfile (str): Path to the command file for sending instructions to MATLAB.
        statusfile (str): Path to the status file indicating MATLAB server readiness.
        process (subprocess.Popen): The MATLAB server process.
        verbose (bool): Flag to enable verbose logging.
        output_queue (queue.Queue): Queue for collecting MATLAB output lines.
        running (bool): Flag to control output reading threads.
        comm (mpi4py.MPI.Comm or None): MPI communicator for rank checking.
    """

    def __init__(self, color=0, matlab_path="matlab", cmdfile="cmdfile", statusfile="statusfile", verbose=False, comm=None):
        """Initialize the MATLAB server configuration.

        Args:
            color (int, optional): Identifier for file naming to support multiple instances. Defaults to 0.
            matlab_path (str, optional): Path to the MATLAB executable. Defaults to "matlab".
            cmdfile (str, optional): Base name for the command file. Defaults to "cmdfile".
            statusfile (str, optional): Base name for the status file. Defaults to "statusfile".
            verbose (bool, optional): Enable verbose logging. Defaults to False.
            comm (mpi4py.MPI.Comm, optional): MPI communicator for rank checking. Defaults to None.
        """
        self.matlab_path = matlab_path
        self.cmdfile = os.path.abspath(f"{cmdfile}_{color}.txt")
        self.statusfile = os.path.abspath(f"{statusfile}_{color}.txt")
        self.process = None
        self.verbose = verbose
        self.output_queue = queue.Queue()  # Queue for asynchronous output handling
        self.running = False  # Controls output reading threads
        self.comm = comm  # MPI communicator

    def kill_matlab_processes(self):
        """Terminate all non-GUI MATLAB processes.

        Scans running processes to identify and terminate non-GUI MATLAB instances based on
        command-line flags. GUI instances are skipped to avoid interfering with user sessions.
        Only rank zero prints output if a communicator is provided.

        Raises:
            Exception: If an error occurs during process termination, it is caught and logged.
        """
        matlab_count = 0
        for proc in psutil.process_iter(['name', 'pid', 'cmdline']):
            try:
                # Check for MATLAB processes (name varies by OS)
                if 'matlab' in proc.info['name'].lower() or 'MATLAB' in proc.info['name']:
                    # if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
                    #     print(f"Found MATLAB process: {proc.info['name']} (PID: {proc.info['pid']})")
                    
                    # Determine if the process is a GUI instance
                    cmdline = proc.info['cmdline']
                    is_gui = True  # Assume GUI unless non-GUI flags are present
                    if cmdline and any(flag in cmdline for flag in ['-nodisplay', '-nodesktop']):
                        is_gui = False
                    
                    if not is_gui:
                        # Terminate non-GUI process
                        if platform.system() == "Windows":
                            proc.terminate()  # Windows uses terminate
                        else:
                            proc.send_signal(signal.SIGTERM)  # Unix uses SIGTERM
                        matlab_count += 1
                        # if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
                        #     print(f"Terminated MATLAB process (PID: {proc.info['pid']})")
                    else:
                        if self.comm is None or self.comm.Get_rank() == 0 and self.verbose:
                            # print(f"Skipped GUI MATLAB process (PID: {proc.info['pid']})")
                            pass
                        
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                continue  # Skip processes that cannot be accessed or are invalid
    
        # if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
        #     if matlab_count == 0:
        #         print("No non-GUI MATLAB processes found to terminate.")
        #     else:
        #         print(f"Terminated {matlab_count} non-GUI MATLAB process(es).")

    def _read_stream(self, stream, stream_name):
        """Read a stream line by line and put lines into the output queue.

        Runs in a separate thread to read MATLAB stdout or stderr in real-time.

        Args:
            stream (file): The stream to read (e.g., process.stdout or process.stderr).
            stream_name (str): Name of the stream ("stdout" or "stderr") for logging.
        """
        while self.running:
            try:
                line = stream.readline().decode('utf-8', errors='ignore').strip()
                if line:
                    self.output_queue.put((stream_name, line))
            except Exception as e:
                self.output_queue.put((stream_name, f"Error reading {stream_name}: {e}"))
                break
            if stream.closed:
                break

    def _print_output(self):
        """Print output from the queue when verbose is enabled.

        Runs in a separate thread to process and print output lines from the queue.
        Only rank zero prints output if a communicator is provided.
        """
        while self.running:
            try:
                stream_name, line = self.output_queue.get(timeout=0.1)
                if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
                    print(f"[ICESEE ⇄ MATLAB {stream_name}] {line}")
                else:
                    from mpi4py import MPI
                    comm = MPI.COMM_WORLD
                    rank = comm.Get_rank()
                    if (self.comm is None or comm.Get_rank() == 0):
                        print(f"[ICESEE ⇄ MATLAB {stream_name}] {line}")
            except queue.Empty:
                continue  # No output available, keep checking

    def launch(self):
        """Launch the MATLAB server and wait for it to be ready.

        Starts the MATLAB server in non-GUI mode, redirects its output, and waits for
        the status file to indicate readiness. Output is printed in real-time if verbose
        and only by rank zero if a communicator is provided.

        Raises:
            SystemExit: If the server fails to start or returns an unexpected status.
            Exception: General errors during launch are caught and logged.
        """
        try:
            self.kill_matlab_processes()
        except Exception as e:
            if self.comm is None or self.comm.Get_rank() == 0:
                print(f"An error occurred: {e}")
            
        # if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
        #     print("[ICESEE ⇄ Launcher] Starting MATLAB server...")
        #     print(f"[ICESEE ⇄ Launcher] Command file: {self.cmdfile}")
        #     print(f"[ICESEE ⇄ Launcher] Status file: {self.statusfile}")
            
        # Clean up old command and status files
        for f in [self.cmdfile, self.statusfile]:
            if os.path.exists(f):
                os.remove(f)
        
        try:
            # Launch MATLAB with non-GUI flags and redirect I/O
            # matlab_cmd = f"{self.matlab_path} -nodesktop -nodisplay -nosplash -nojvm -r \"matlab_server('{self.cmdfile}', '{self.statusfile}')\""
            matlab_cmd = f"{self.matlab_path} -nodesktop -nosplash -r \"matlab_server('{self.cmdfile}', '{self.statusfile}')\""

            self.process = subprocess.Popen(
                matlab_cmd,
                shell=True,
                stdout=subprocess.PIPE,  # Redirect stdout
                stderr=subprocess.PIPE,  # Redirect stderr
                stdin=subprocess.PIPE,   # Redirect stdin
                preexec_fn=os.setsid    # Create new process group for signal handling
            )
            
            self.running = True

            # Start threads to handle real-time output
            stdout_thread = threading.Thread(target=self._read_stream, args=(self.process.stdout, "stdout"))
            stderr_thread = threading.Thread(target=self._read_stream, args=(self.process.stderr, "stderr"))
            output_thread = threading.Thread(target=self._print_output)

            stdout_thread.daemon = True
            stderr_thread.daemon = True
            output_thread.daemon = True

            stdout_thread.start()
            stderr_thread.start()
            output_thread.start()
            
            # Wait for server to signal readiness via status file
            timeout = 1000  # seconds
            start_time = time.time()
            while not os.path.exists(self.statusfile):
                if time.time() - start_time > timeout:
                    if self.comm is None or self.comm.Get_rank() == 0:
                        print("[ICESEE ⇄ Launcher] Error: MATLAB server failed to start within timeout.")
                    self.running = False
                    os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                    sys.exit(1)
                time.sleep(0.5)
            
            # Check server status
            with open(self.statusfile, 'r') as f:
                status = f.read().strip()
            if status != 'ready':
                if self.comm is None or self.comm.Get_rank() == 0:
                    print(f"[ICESEE ⇄ Launcher] Error: Unexpected status '{status}'.")
                self.running = False
                os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                sys.exit(1)
            
            if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
                print("[ICESEE ⇄ Launcher] MATLAB server is ready.")
        except Exception as e:
            if self.comm is None or self.comm.Get_rank() == 0:
                print(f"[ICESEE ⇄ Launcher] Error launching MATLAB server: {e}")
            self.running = False
            if self.process:
                os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
            sys.exit(1)

    # def send_command(self, command, timeout=3600):
    #     """Send a command to MATLAB and wait for it to be processed.

    #     Writes the command to the command file and waits for MATLAB to process it
    #     (indicated by the file being deleted). Provides periodic status updates if verbose
    #     and only by rank zero if a communicator is provided. The verbose interval scales
    #     dynamically with the timeout.

    #     Args:
    #         command (str): The MATLAB command to execute.
    #         timeout (int, optional): Maximum time to wait for command processing (seconds). Defaults to 3600.

    #     Returns:
    #         bool: True if the command was processed successfully, False if it timed out.
    #     """
    #     # Set dynamic verbose interval: at least 60 seconds, or timeout / 20
    #     verbose_interval = max(60.0, timeout / 20.0)

    #     if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
    #         print(f"[ICESEE ⇄ Launcher] Sending command: {command}")
        
    #     try:
    #         with open(self.cmdfile, 'w') as f:
    #             f.write(command)
    #     except OSError as e:
    #         if self.comm is None or self.comm.Get_rank() == 0:
    #             print(f"[ICESEE ⇄ Launcher] Error: Failed to write command file: {e}")
    #         return False
        
    #     # Wait for command to be processed (file deleted)
    #     start_time = time.time()
    #     sleep_time = 0.2  # Initial sleep interval
    #     max_sleep = 60.0  # Maximum sleep interval
    #     last_verbose_time = start_time  # Track last verbose print
    #     warning_threshold = timeout * 0.8  # Warn at 80% of timeout
    #     warning_issued = False
        
    #     while os.path.exists(self.cmdfile):
    #         elapsed_time = time.time() - start_time
    #         if elapsed_time > timeout:
    #             if self.comm is None or self.comm.Get_rank() == 0:
    #                 print(f"[ICESEE ⇄ Launcher] Error: Command execution timed out after {timeout} seconds.")
    #             return False
            
    #         # Issue warning if approaching timeout
    #         if not warning_issued and elapsed_time > warning_threshold:
    #             if self.comm is None or self.comm.Get_rank() == 0:
    #                 print(f"[ICESEE ⇄ Launcher] Warning: Command has been running for {elapsed_time:.1f}s, approaching timeout of {timeout}s.")
    #             warning_issued = True
            
    #         # Print periodic status if verbose
    #         if self.verbose and (self.comm is None or self.comm.Get_rank() == 0) and (time.time() - last_verbose_time) >= verbose_interval:
    #             print(f"[ICESEE ⇄ Launcher] Waiting for command to be processed... ({elapsed_time:.1f}s elapsed)")
    #             last_verbose_time = time.time()
            
    #         time.sleep(sleep_time)
    #         # Gradually increase sleep time to reduce CPU usage
    #         sleep_time = min(sleep_time + (elapsed_time / 20.0), max_sleep)
    #         if sleep_time == max_sleep and not warning_issued:
    #             if self.comm is None or self.comm.Get_rank() == 0:
    #                 print("[ICESEE ⇄ Launcher] Warning: Slow command processing detected.")
        
    #     if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
    #         print("[ICESEE ⇄ Launcher] Command processed successfully.")
    #     return True

    def send_command(self, command, timeout=3600):
        """Send a command to MATLAB and wait for it to be processed.

        Writes the command to the command file and waits for MATLAB to process it
        (indicated by the file being deleted). Provides periodic status updates if verbose
        and only by rank zero if a communicator is provided. The verbose interval scales
        dynamically with the timeout.

        Args:
            command (str): The MATLAB command to execute.
            timeout (int, optional): Maximum time to wait for command processing (seconds). Defaults to 3600.

        Returns:
            bool: True if the command was processed successfully, False if it timed out.
        """
        # Set dynamic verbose interval: at least 60 seconds, or timeout / 20
        verbose_interval = max(60.0, timeout / 20.0)

        if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
            print(f"[ICESEE ⇄ Launcher] Sending command: {command}")
        
        try:
            with open(self.cmdfile, 'w') as f:
                f.write(command)
        except OSError as e:
            if self.comm is None or self.comm.Get_rank() == 0:
                print(f"[ICESEE ⇄ Launcher] Error: Failed to write command file: {e}")
            return False
        
        # Wait for command to be processed (file deleted)
        start_time = time.time()
        sleep_time = 0.2  # Initial sleep interval
        min_sleep = 0.1   # Minimum sleep interval to prevent busy-waiting
        max_sleep = 60.0  # Maximum sleep interval
        last_verbose_time = start_time  # Track last verbose print
        warning_threshold = timeout * 0.8  # Warn at 80% of timeout
        warning_issued = False
        last_loop_start = start_time  # Track start of the current loop iteration
        
        while os.path.exists(self.cmdfile):
            elapsed_time = time.time() - start_time
            if elapsed_time > timeout:
                if self.comm is None or self.comm.Get_rank() == 0:
                    print(f"[ICESEE ⇄ Launcher] Error: Command execution timed out after {timeout} seconds.")
                return False
            
            # Issue warning if approaching timeout
            if not warning_issued and elapsed_time > warning_threshold:
                if self.comm is None or self.comm.Get_rank() == 0:
                    print(f"[ICESEE ⇄ Launcher] Warning: Command has been running for {elapsed_time:.1f}s, approaching timeout of {timeout}s.")
                warning_issued = True
            
            # Print periodic status if verbose
            if self.verbose and (self.comm is None or self.comm.Get_rank() == 0) and (time.time() - last_verbose_time) >= verbose_interval:
                print(f"[ICESEE ⇄ Launcher] Waiting for command to be processed... ({elapsed_time:.1f}s elapsed)")
                last_verbose_time = time.time()
            
            # Measure the time taken by the previous loop iteration
            current_time = time.time()
            loop_duration = current_time - last_loop_start
            last_loop_start = current_time  # Update for the next iteration
            
            # Set sleep_time to the previous loop's duration, bounded by min_sleep and max_sleep
            sleep_time = min(max(loop_duration, min_sleep), max_sleep)
            
            time.sleep(sleep_time)
            
            # Issue warning if sleep_time reaches max_sleep
            if sleep_time == max_sleep and not warning_issued:
                if self.comm is None or self.comm.Get_rank() == 0:
                    print("[ICESEE ⇄ Launcher] Warning: Slow command processing detected.")
        
        if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
            print("[ICESEE ⇄ Launcher] Command processed successfully.")
        return True


    def shutdown(self):
        """Attempt to gracefully shut down the MATLAB server.

        Sends an 'exit' command to MATLAB and waits for termination. Captures any
        remaining output and handles forced termination if necessary. Only rank zero
        prints output if a communicator is provided.

        Raises:
            subprocess.TimeoutExpired: If the process does not terminate within the timeout.
        """
        if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
            print("[ICESEE ⇄ Launcher] Attempting to shut down MATLAB server...")

        self.running = False  # Stop output threads

        if self.send_command("exit"):
            try:
                # Wait for process to terminate and collect remaining output
                stdout, stderr = self.process.communicate(timeout=5)
                if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
                    if stdout:
                        print("[ICESEE⇄MATLAB]", stdout.decode('utf-8', errors='ignore'))
                    if stderr:
                        print("[ICESEE⇄MATLAB]", stderr.decode('utf-8', errors='ignore'))
                    print("[ICESEE ⇄ Launcher] MATLAB server shut down successfully.")
            except subprocess.TimeoutExpired:
                if self.comm is None or self.comm.Get_rank() == 0:
                    print("[ICESEE ⇄ Launcher] Warning: MATLAB process did not terminate in time, forcing termination.")
                os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                self.process.wait(timeout=300)
                if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
                    print("[ICESEE ⇄ Launcher] MATLAB server terminated.")
        else:
            if self.comm is None or self.comm.Get_rank() == 0:
                print("[ICESEE ⇄ Launcher] Error: Failed to shut down MATLAB server gracefully.")
            os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
            self.process.wait(timeout=5)
            if self.comm is None or self.comm.Get_rank() == 0:
                print("[ICESEE ⇄ Launcher] MATLAB server terminated.")

    def reset_terminal(self):
        """Reset terminal settings to restore normal behavior, if applicable.

        Skips reset in non-interactive or MPI environments to avoid interference.
        Uses 'stty sane' to restore terminal settings. Only rank zero prints output
        if a communicator is provided.

        Raises:
            subprocess.CalledProcessError: If the terminal reset command fails.
        """
        if not sys.stdin.isatty() or any(key in os.environ for key in ["MPIEXEC", "OMPI_COMM_WORLD_RANK", "PMI_RANK", "SLURM_JOB_ID"]):
            if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
                print("[ICESEE ⇄ Launcher] Skipping terminal reset (non-interactive or MPI environment).", file=sys.stderr, flush=True)
            return
        try:
            subprocess.run(['stty', 'sane'], check=True)
            if self.verbose and (self.comm is None or self.comm.Get_rank() == 0):
                print("[ICESEE ⇄ Launcher] Terminal settings reset successfully.", file=sys.stderr, flush=True)
        except subprocess.CalledProcessError as e:
            if self.comm is None or self.comm.Get_rank() == 0:
                print(f"[ICESEE ⇄ Launcher] Warning: Failed to reset terminal settings: {e}", file=sys.stderr, flush=True)

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
                    verbose=self.verbose,
                    comm=self.comm
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
            print("------ ICESEE⇄MATLAB ------")
            print(trimmed_stdout.strip())

            if stderr.strip():
                print("------ ICESEE⇄MATLAB ------")
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
#             print("------ ICESEE⇄MATLAB ------")
#             print(trimmed_stdout.strip())

#             if stderr.strip():  # Only print stderr if there's content
#                 print("------ ICESEE⇄MATLAB ------")
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

    # print(f"[ICESEE] Added ISSM directory and subdirectories from path: {issm_dir}")


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


def setup_ensemble_data(Nens, base_data_dir='./Models/ens_id_0', base_kwargs_file='model_kwargs_0.mat'):
    import os
    import shutil
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    base_data_dir = os.path.abspath(base_data_dir)
    base_kwargs_file = os.path.abspath(base_kwargs_file)
    
    if rank == 0:
        if not os.path.isdir(base_data_dir):
            raise FileNotFoundError(f"[Rank {rank}] Base directory {base_data_dir} not found")
        if not os.path.isfile(base_kwargs_file):
            raise FileNotFoundError(f"[Rank {rank}] Base kwargs file {base_kwargs_file} not found")
            
        for ens in range(Nens):
            ens_dir = os.path.abspath(f'./Models/ens_id_{ens}')
            kwargs_file = f'model_kwargs_{ens}.mat'
            
            if ens != 0:  # Skip ens_id_0 (base directory)
                if os.path.exists(ens_dir):
                    shutil.rmtree(ens_dir)
                os.makedirs(ens_dir)
                
                for root, _, files in os.walk(base_data_dir):
                    rel_path = os.path.relpath(root, base_data_dir)
                    os.makedirs(os.path.join(ens_dir, rel_path), exist_ok=True)
                    for file_name in files:
                        os.link(os.path.join(root, file_name), os.path.join(ens_dir, rel_path, file_name))
                
                if os.path.exists(kwargs_file):
                    os.remove(kwargs_file)
                os.link(base_kwargs_file, kwargs_file)
    
    comm.Barrier()
    
    ens_id = rank
    if ens_id < Nens:
        ensemble_dir = os.path.abspath(f'./Models/ens_id_{ens_id}')
        ensemble_kwargs = f'model_kwargs_{ens_id}.mat'
        if not os.path.isdir(ensemble_dir) or not os.path.isfile(ensemble_kwargs):
            raise FileNotFoundError(f"[Rank {rank}] Cannot access {ensemble_dir} or {ensemble_kwargs}")
        return ensemble_dir, ensemble_kwargs
    return None, None


def setup_reference_data(reference_data_dir, reference_data, use_reference_data=True):
    """
    Create ensemble directories with hard-linked reference data file for read-only access.
    
    Parameters:
    - reference_data_dir: Directory containing the reference data file.
    - reference_data: Name of the reference data file.
    - use_reference_data: Flag to enable/disable reference data setup.
    
    Returns:
    - rank_data_dir: Path to the rank's ensemble directory (e.g., './Models/ens_id_X').
    - rank_data_file: Path to the rank's reference data file.
    """
    from mpi4py import MPI
    import os
    import shutil

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    initial_data = os.path.abspath(os.path.join(reference_data_dir, reference_data))
    rank_data_dir = os.path.abspath(f'./Models/ens_id_{rank}')
    rank_data_file = os.path.join(rank_data_dir, reference_data)

    if use_reference_data and rank == 0:
        # Verify reference data exists
        if not os.path.isfile(initial_data):
            raise FileNotFoundError(f"[Rank {rank}] Reference data {initial_data} not found")

        # Create directories and hard link the reference file
        for _rank in range(size):
            rank_data_dir = os.path.abspath(f'./Models/ens_id_{_rank}')
            link_path = os.path.join(rank_data_dir, reference_data)

            # Create or clear directory
            if os.path.exists(rank_data_dir):
                if os.path.isdir(rank_data_dir):
                    shutil.rmtree(rank_data_dir)
                else:
                    os.remove(rank_data_dir)
            os.makedirs(rank_data_dir, exist_ok=True)

            # Hard link the reference file
            if os.path.exists(link_path):
                os.remove(link_path)
            try:
                os.link(initial_data, link_path)
            except OSError as e:
                raise RuntimeError(f"[Rank {rank}] Failed to create hard link {link_path} -> {initial_data}: {e}")

        # print(f"[Rank 0] Created {size} ensemble directories with hard-linked reference data")

    # Synchronize all ranks
    comm.Barrier()

    # Verify access for this rank
    if use_reference_data:
        if not os.path.isdir(rank_data_dir) or not os.path.isfile(rank_data_file):
            raise FileNotFoundError(f"[Rank {rank}] Cannot access {rank_data_dir} or {rank_data_file}")
        return rank_data_dir, rank_data_file
    return None, None


# make data available in all ensemble directories
def setup_ensemble_intial_data(Nens, reference_data_dir, reference_data):
    """
    Create ensemble directories with hard-linked reference data file for read-only access.
    
    Parameters:
    - reference_data_dir: Directory containing the reference data file.
    - reference_data: Name of the reference data file.
    - use_reference_data: Flag to enable/disable reference data setup.
    
    Returns:
    - rank_data_dir: Path to the rank's ensemble directory (e.g., './Models/ens_id_X').
    - rank_data_file: Path to the rank's reference data file.
    """
    from mpi4py import MPI
    import os
    import shutil

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    initial_data = os.path.abspath(os.path.join(reference_data_dir, reference_data))
    rank_data_dir = os.path.abspath(f'./Models/ens_id_{rank}')
    rank_data_file = os.path.join(rank_data_dir, reference_data)

    if rank == 0:
        # Verify reference data exists
        if not os.path.isfile(initial_data):
            raise FileNotFoundError(f"[Rank {rank}] Reference data {initial_data} not found")

        # Create directories and hard link the reference file
        for ens in range(Nens):

            # skip ens_id_0 (base directory)
            if ens == 0:
                continue

            rank_data_dir = os.path.abspath(f'./Models/ens_id_{ens}')
            link_path = os.path.join(rank_data_dir, reference_data)

            # Hard link the reference file
            if os.path.exists(link_path):
                os.remove(link_path)
            try:
                os.link(initial_data, link_path)
            except OSError as e:
                raise RuntimeError(f"[Rank {rank}] Failed to create hard link {link_path} -> {initial_data}: {e}")





           