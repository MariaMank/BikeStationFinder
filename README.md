# BikeStationFinder

The application allows you to find the best location for new bike-sharing stations using Location-Allocation analysis.

In the BaseClass.py file you can find the base class which allows you to import OSM data, save it chosen format and join with different GeoDataFrame. 
In preprossecing.py file you can find classes and functions responsible for preparation of all data needed for Network Analysis (Geocoding of the area of interest, creating propositions of the new location of bike sharing stations, creating the Footway Net, Choosing and saving all Points Of Interest relevant to the study subject)
In solver.py file you can find implementation of the class which responsible for the actual spatial analysis (Location - Allocation analysis using PySal library).
In the main file.py you can find example of calling the solver of BikeStationFinder.
![RozneWagiUrsus700](https://github.com/user-attachments/assets/a39d5f81-135c-47ba-b735-88da35052383)





Note: The application uses Polish Boundary System WFS - will not work for the area outside of Poland without adjustments
