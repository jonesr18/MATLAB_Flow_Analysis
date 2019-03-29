# MATLAB_Flow_Analysis
This repository is for object-oriented MATLAB-based flow cytometry data analysis. We provide a number of static classes and objects for use in processing and analyzing multi-color samples. 

Static classes act like namespaces; to use them, write out something like the following:

    plt = Plotting();
    plt.binHeatmap( ... );

The repository also contains a data structure object called FlowData(), which is built around a struct-based storage of flow data with the purpose being easier processing and data management, specifically with the following goals:

1. Connect sample IDs to actual experimental parameters
2. Trivialize combining gates/different data types
3. Minimize memory overhead

Built-in capabilities include:
* Gating (thresholds, polygons)
* Matrix-based color compensation & autofluorescence subtraction
* Conversion to MEFL units
* Biexponential transformations
* Many plotting tools, including biexponential axes scaling/labeling
* Data modeling/fitting tools
* Fluorescence binning
* Calculation of sample statistics
