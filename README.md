# Matrix-Transposition-Optimization

## Instructions for Running on Windows

### Step 1: Install MobaXterm

1. Download **MobaXterm Home Edition v24.3**:
   - Go to [MobaXterm Official Website](https://mobaxterm.mobatek.net/).
   - Click on **GET MOBAXTERM NOW!**.
   - Select the **Free Edition** and download the **MobaXterm Home Edition v24.3 (Installer edition)**.

2. Install MobaXterm:
   - Unzip the downloaded folder.
   - Open the unzipped folder and double-click on the file: `MobaXterm_installer_24.3`.
   - Follow the installation steps:
     - Click **Next** on the first screen.
     - Accept the terms and click **Next**.
     - Keep the default folder and click **Next** again.
     - Click **Install**.
   - Grant permission to modify your computer by clicking **Yes**.
   - When the installation finishes, click **Finish**.

   > **Note:** We will open MobaXterm later in the process.

---

### Step 2: Install GlobalProtect VPN

1. Download **GlobalProtect**:
   - Visit the [GlobalProtect download page](https://vpn.youruniversitydomain.com/).
   - Select the version compatible with your system (32-bit or 64-bit) based on the instructions on the page.

2. Install GlobalProtect:
   - Run the downloaded file.
   - Click **Next** on the three setup screens to proceed.
   - Grant permission to modify your computer by clicking **Yes**.
   - Once the installation is complete, click **Close**.

---

### Step 3: Connect to the University Network via VPN

1. Locate the **GlobalProtect** icon:
   - Click the up arrow on the bottom-right corner of the taskbar.
   - Look for the **grey globe** icon and click on it.

2. Connect to the VPN:
   - In the **Portal** field, type or paste: `vpn-mfa.icts.unitn.it`.
   - Click **Connect**.
   - When prompted, log in with your university account credentials or use **SPID** if requested.

---

### Step 4: Proceed with GitHub and HPC Setup

When the VPN is activated, proceed with the following steps:

1. Open the GitHub link: [Matrix Transposition Optimization](https://github.com/ale-bena/Matrix-Transposition-Optimization.git).
2. Click the green **Code** button and select **Download ZIP**.
3. Unzip the folder.

---

### Step 5: Configure MobaXterm for HPC Access

1. **Open MobaXterm** and start a new local session.
2. If Windows Defender prompts you to choose network settings, select **Private Networks** for safety.

3. Connect to the HPC cluster by typing the command:
   ```bash
   ssh name.surname@hpc.unitn.it
 Replace `name` and `surname` with your personal credentials.
   - Enter your account password when prompted.
   - It is recommended not to save the password for security reasons—close the pop-up window instead.

2. **Navigate Your HPC Directory**
   - Once logged in, check the left panel in **MobaXterm** to see your directory files.
   - If you don’t see the directory, click the **orange world icon** to display it.

3. **Modify the `.bashrc` File**
   - Locate the `.bashrc` file in your directory.
   - Open the file:
     - Double-click it, or
     - Right-click and select **Open with default text editor**.
   - Add the following lines under the section `# User specific aliases and functions`:
     ```bash
     module load gcc91
     alias gcc=gcc-9.1.0
     ```

4. **Copy Project Files to HPC**
   - Download the project folder from GitHub.
   - Drag and drop the folder into the **MobaXterm** file browser. The files will be automatically copied to your HPC environment.

5. **Edit and Run the PBS File**
   - Navigate to the project directory:
     ```bash
     cd project
     ```
   - Open the `exec.pbs` file:
     - Right-click the file and select **Open with default text editor**.
     - Locate the line:
       ```bash
       # Select the working directory (change this line)
       cd /home/alessandro.benassi/project
       ```
     - Replace `alessandro.benassi` with your `name.surname`.

   - Run the PBS file to execute the project:
     ```bash
     qsub exec.pbs
     ```
After some time, a folder named `result_1` should appear, and you should see the files `times_tables.csv` and `times_average.csv` in the project folder. When the execution terminates, these files will be moved to the `result_1` directory.

6. **Project Structure Overview**
   - The project contains multiple implementations:
     - **Sequential**: Three files (header, function, and main).
     - **Implicit**: Three files (header, function, and main).
     - **OpenMP**: Three files (header, function, and main).
     - **Block-based OpenMP**: Three files (header, function, and main).
   - This modular approach ensures reusability and clarity for external users.

7. **Compute Speedup**
   - Run the **Sequential** implementation as a baseline.
   - Use the results to compute speedup for the **Implicit** and **OpenMP** implementations.
   - Without sequential execution data, speedup cannot be calculated.

---

Follow these instructions carefully to complete the setup and run the project on the HPC cluster.

Note: the comments of the functions are presents in the sequential code. The other files have the comments only for the functions that are different.

