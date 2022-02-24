# Standard_vs_KDE_histogram
Traditional histograms vs Kernel Density Estimators (KDEs).

## **What is this repository for?**

In this script, I compare the properties of traditional histograms to those from KDEs. I present my own function to build a KDE-histogram with a 1D gaussian as kernel and show how these techniques can be used to draw information from data. Besides the different techniques, I also record some other plotting routines e.g. scatter, plot, spine, color bars, etc. The code was developed to analyze spectroscopic data e.g. spectroscopic redshifts from MUSE data. 

## **Installing QtCounterpart**

No installation is needed. Simply cloning the GitHub repository and importing the script is enough. Hence, 

```
    git clone https://github.com/YohanaHerrero/Standard_vs_KDE_histogram.git
```

The code is written in Python and uses a range of default packages included in standard installations of Python:

### **Standard packages:**

- numpy  
- matplotlib
- math
- seaborn

### **Special packages:**

- astropy 

After adding the Standard_vs_KDE_histogram directory to the PYTHONPATH or changing location to the Standard_vs_KDE_histogram directory, the repository can be imported in python with:

```
    import Standard_vs_KDE_histogram
```

Besides the python packages, you will also need a fits table containing the following columns:

- RA in degrees
- DEC in degrees
- Z 

Instead, the three parameters above can be replaced by various kind of data, according to specific needs.
