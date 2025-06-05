# Software Installation Guide

If you're using a work laptop, please get permission from your local IT helpdesk to get admin access or consider using a personal computer. Once you have completed the installation, don't forget to update the Google form if you have previously filled it out.

## 1. AWS CLI (Command Line Interface)

   - **For Windows**:
     1. Download the installer from [AWS CLI website](https://aws.amazon.com/cli/).
     2. Run the installer and follow the prompts.
     3. Verify installation by opening `Command Prompt` and typing:
        ```bash
        aws --version
        ```

   - **For macOS**:
     1. Open Terminal and use Homebrew to install:
        ```bash
        brew install awscli
        ```
     2. Verify installation:
        ```bash
        aws --version
        ```

   - **For Linux**:
     1. Use the following commands:
        ```bash
        curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
        unzip awscliv2.zip
        sudo ./aws/install
        ```
     2. Verify installation:
        ```bash
        aws --version
        ```

## 2. R (Statistical Programming Language)

   - **For Windows**:
     1. Download R from the [official website](https://cran.r-project.org/).
     2. Run the installer and follow the prompts.

   - **For macOS**:
     1. Download R from [CRAN for macOS](https://cran.r-project.org/bin/macosx/).
     2. Open the `.pkg` file and follow the instructions.

   - **For Linux**:
     1. Use the following command (Ubuntu example):
        ```bash
        sudo apt install r-base
        ```
     2. Verify installation by typing:
        ```bash
        R
        ```

## 3. FastK v1.1

   - **For all platforms**:
     1. Download the binary from the [FastK GitHub release page](https://github.com/thegenemyers/FASTK).
     2. Unzip the file and place it in a directory of your choice.
     3. Add that directory to your `PATH` or run it directly from the terminal.

## 4. Genomescope v2.0.1

   - **For all platforms**:
     1. Go to the [Genomescope repository](https://github.com/schatzlab/genomescope).
     2. Download the appropriate version or install via `conda`:
        ```bash
        conda install -c bioconda genomescope
        ```
     3. Verify installation by running:
        ```bash
        genomescope
        ```

## 5. sourmash v4.9.0
> If you don't have `conda` or `mamba` installed, please set up `conda` as described below.

   - **For all platforms**:
     1. Install via `conda`:
        ```bash
         conda install -c conda-forge -c bioconda sourmash
        ```
     2. Verify installation by running:
        ```bash
        sourmash --version
        ```

# Setting Up Conda for Bioinformatics Software

### Step 1: Install Miniforge on Windows or macOS

**Miniforge** is a lightweight Conda installer that defaults to using `conda-forge`, which is free and community-maintained. This avoids commercial restrictions tied to the `defaults` channel (e.g., Anaconda).

👉 Download Miniforge from the official releases page:  
[Miniforge Releases](https://github.com/conda-forge/miniforge/releases/latest)

Choose the appropriate installer:
- **Windows**: `Miniforge3-Windows-x86_64.exe`
- **macOS (Intel)**: `Miniforge3-MacOSX-x86_64.sh`
- **macOS (M1/M2 ARM)**: `Miniforge3-MacOSX-arm64.sh`

### 🖥️ Step 2: Install Miniforge (Example for Linux)

**For Linux users only**:

Install `wget` if not available:
```bash
sudo apt update && sudo apt install -y wget
```

Download the installer:

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
```

Run the installer:

```bash
bash Miniforge3-Linux-x86_64.sh
```

Follow the prompts:

Press ENTER to begin

Accept the license

Choose install path (default: ~/miniforge3)

Let it initialize Conda when prompted

### Step 3: Initialize and Verify Conda
After installation, either:

```bash
source ~/miniforge3/bin/activate
```

Or restart your terminal.

Check Conda version:

```bash
conda --version
```

Expected output:

```bash
conda 24.11.3
```

🛠️ Additional (Optional): Manually Add Conda to Your Shell
If Conda isn’t working after install, run:

```bash
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
```

For Zsh users:

```bash
eval "$($HOME/miniforge3/bin/conda shell.zsh hook)"
```

Then initialize Conda for future sessions:

```bash
conda init
```
This updates your ~/.bashrc or ~/.zshrc file to load Conda automatically.

Restart the terminal, or:

```bash
source ~/.bashrc    # for bash
source ~/.zshrc     # for zsh
```

### Step 4: Configure Conda Channels
By default, Miniforge uses conda-forge, but for bioinformatics you’ll want to add bioconda.

Check current channels:

```bash
conda config --show channels
```

Remove paid/proprietary channels:

```bash
conda config --remove channels defaults
conda config --remove channels anaconda
conda config --remove channels r
conda config --remove channels pro
```

Add trusted community channels:

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

Confirm configuration:

```bash
conda config --show channels
```

You should see:

channels:
  - conda-forge
  - bioconda
🔧 Alternative: Edit .condarc manually

```bash
nano ~/.condarc
```

Paste this:

```bash
channels:
  - conda-forge
  - bioconda
  - nodefaults
denylist_channels: [defaults, anaconda, r, main, pro] #!final
```

### Step 5: Install softwares

```bash
conda install -c bioconda fastk
conda install -c bioconda genomescope
```

---

If you have any questions or encounter issues during installation, feel free to ask ChatGPT or reach out via the Question and Comments on Learning and Management System (LMS)

---
