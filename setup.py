import subprocess
import sys

# List of required dependencies
dependencies = [
    'csv',
    'time',
    'pandas',
    'scipy',
    'scikit-image',
    'imageio',
    'numpy',
    'matplotlib',
    'matplotlib.mlab',
    'nd2',
    'warnings'
    'openpyxl',
]

def install_dependencies():
    for package in dependencies:
        try:
            # Check if the package is already installed
            __import__(package)
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])

if __name__ == "__main__":
    print("Checking and installing dependencies...")
    install_dependencies()
    print("Dependencies installed successfully.")
