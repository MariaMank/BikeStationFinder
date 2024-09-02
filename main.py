from solver import BikeStationSolver

if __name__ == '__main__':
    algorytmy = ['MCLP', 'LSCP', 'P-median', 'P-center', 'All']
    results_path = f"./dane/"
    propozycje_path = results_path + 'propozycjeP.geojson'
    istniejace_path = results_path + 'existing_stations.geojson'
    siec = results_path + 'NetP.geojson'
    demands = "./mclpF/UrsusDane/living_houses.geojson"



    solver = BikeStationSolver(500,
                                3,
                                800,
                                area='Brzeźno',
                                algorytm='MCLP')
    solver.solve()
    solver.getResults(results_path)
    """
    solver = BikeStationSolver(500,
                               3,
                               200,
                               'All',
                               existing=istniejace_path,
                               proposiotions=propozycje_path,
                               net=siec,
                               pois=demands,
                               weights='wagi')

    """
