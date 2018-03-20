# MATLAB_Flow_Analysis
MATLAB Flow Analysis Code from the Weiss Lab at MIT

This code compiles many useful functions from various lab members into Static Classes which act like namespaces. To use them, write out something like the following:

    plt = Plotting();
    plt.binHeatmap( ... );

It also contains a data structure called FlowData, which is built around a struct-based storage of flow data with the purpose being easier processing and data management, specifically with the following goals:

1. Connect sample IDs to actual experimental parameters
2. Trivialize combining gates/different data types
3. Minimize memory overhead
