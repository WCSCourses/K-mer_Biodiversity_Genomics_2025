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
     1. Download the binary from the [FastK GitHub release page](https://github.com/axel4/FastK).
     2. Unzip the file and place it in a directory of your choice.
     3. Add that directory to your `PATH` or run it directly from the terminal.

## 4. Genomescope v2.0.1

   - **For all platforms**:
     1. Go to the [Genomescope repository](https://github.com/bcgsc/Genomescope).
     2. Download the appropriate version or install via `conda`:
        ```bash
        conda install -c bioconda genomescope
        ```
     3. Verify installation by running:
        ```bash
        genomescope
        ```

---

If you have any questions or encounter issues during installation, feel free to ask ChatGPT or reach out via the Question and Comments on Learning and Management System (LMS)

---
