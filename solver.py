import datetime
import geopandas as gpd
import shapely
import pulp
import warnings
from datetime import datetime
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # ignore deprecation warning - GH pysal/spaghetti#649
    import spaghetti
    import spopt
    from spopt.locate import MCLP
import numpy as np
import pandas as pd
from baseClass import BaseClass
from preprocessing import Net, Pois, Proposition, Bariers, Obszar, Building




class BikeStationSolver:
    def __init__(self, service_radius, solutionsNumber, frequencyPropPoints,algorytm ,area=None, existing=None, net=None,
                 proposiotions=None, weights=None, pois=None):
        self.obszar = Obszar(area)
        self.rad = service_radius
        self.skok = frequencyPropPoints
        self.algorytm = algorytm
        print(datetime.now(), "start")
        print(datetime.now(), "tworze siec")
        self.net = Net(self.obszar, ready=net)
        print(datetime.now(), 'pobieram dane o istniejacych stacjach')
        self.existingStations = BaseClass(self.obszar.area_id, "node[amenity = bicycle_rental](area.b)" ,ready=existing)
        self.existingStations = self.existingStations.dataInFrame
        if self.existingStations.size != 0:
            self.existingStations = self.existingStations.assign(
            predefined_loc=np.ones(len(self.existingStations.geometry), dtype=bool))
            self.solutionNumber = solutionsNumber + len(self.existingStations.geometry)
        print(datetime.now(), 'tworze punkty na potencjalne stacje')
        #else:
        self.bariery = Bariers(self.obszar)  # bariery - tam gdzie nie moze byc stacji - budynki i drogi, tory
        self.propozycje = Proposition(self.skok, self.bariery,
                                          self.net.NetDataFrame, ready=proposiotions)  # kwadratowa siatka propozycji (punkty) po selekcji tego co juz nie moze byc - jakieś 1.5 min
        self.propozycje = self.propozycje.dataInFrame
        self.propozycje = self.propozycje.assign(predefined_loc=np.zeros(len(self.propozycje.geometry), dtype=bool))
        self.propozycje = pd.concat([self.propozycje, self.existingStations], join='inner')

        print(datetime.now(), 'pobieram i waguje budynki mieszkaniowe i generatory ruchu')
        if pois == None:
            self.budynki = Building(self.obszar)
            self.pois = Pois(self.obszar, self.budynki)
            self.pois = self.pois.pois_inDataFrame
            self.budynki = self.budynki.living_houses
            self.weights = self.pois['weights']
        else:
            self.pois = gpd.read_file(pois)
            self.weights = self.pois[weights]
        print(datetime.now(), 'tworze polaczenia propozycji, pois i budynków mieszkaniowych z reszta sieci')
        self.net.createConnectionToTheStations(self.propozycje, self.existingStations, self.pois)
        print('Ilość generatorów:', len(self.pois))
        print('Ilość elementów w sieci:', len(self.net.NetDataFrame))
        print('Suma długości sieci:', np.sum(self.net.NetDataFrame.length))
        print('Ilość propozycji stacji:', len(self.propozycje))
        print('Promień:', self.rad)
        print('Skok:', self.skok)
        print('Obszar:', self.obszar.name)
        print('Algorytm:', self.algorytm)

    def solve(self):
        print(datetime.now(), 'szukam rozwiazania')
        # set the solver
        solver = pulp.COIN_CMD(msg=False, warmStart=True, keepFiles=False)
        #czasem może być konieczne znalezienie plików .exe COIN CMD i dodanie jako argument do powyzszej funkcji
        # path=r'D.\env\Lib\site-packages\pulp\solverdir\cbc\win\64\cbc.exe'
        print(datetime.now(), 'tworze macierz odleglosci po sieci')
        self.cosMatrix = self.createCostsMatrix()
        predef = self.propozycje['predefined_loc'].array
        if self.algorytm == 'MCLP':
            mclp_from_cm = MCLP.from_cost_matrix(
                self.cosMatrix,
                self.weights,
                self.rad,
                p_facilities=self.solutionNumber,
                name=f"mclp-BikeStations-{self.obszar.name}",
                predefined_facilities_arr=predef
            )
            print(datetime.now(), 'zaczynam rozwiazanie')
            self.problem_solution = mclp_from_cm.solve(solver)
        elif self.algorytm == 'LSCP':
            lscp_from_cm = spopt.locate.LSCP.from_cost_matrix(
                cost_matrix=self.cosMatrix,
                service_radius=self.rad,
                name=f"lscp-BikeStations-{self.obszar.name}",
                predefined_facilities_arr=predef
            )
            try:
                self.lscp_from_cm = lscp_from_cm.solve(solver)
            except RuntimeError:
                self.problem_solution = False
                print('Error: Location Set Covering Problem nie ma rozwiazania - sprawdź poprawność wprowadzonych danych')
        elif self.algorytm == 'LSCPB':
            try:
                self.lscpb_from_cm = spopt.locate.LSCPB.from_cost_matrix(
                cost_matrix=self.cosMatrix,
                service_radius=self.rad,
                solver=solver,
                name=f"lscp-BikeStations-{self.obszar.name}",
                predefined_facilities_arr=predef
            )
            except RuntimeError:
                self.problem_solution = False
                print( 'Error: Location Set Covering Problem with Backup nie ma rozwiazania - sprawdź poprawność wprowadzonych danych')
        elif self.algorytm == 'P-median':
            pmed_from_cm = spopt.locate.PMedian.from_cost_matrix(
                cost_matrix=self.cosMatrix,
                weights= self.weights,
                name=f"p-median-BikeStations-{self.obszar.name}",
                p_facilities=self.solutionNumber,
                predefined_facilities_arr=predef
            )
            try:
                self.problem_solution = pmed_from_cm.solve(solver)
            except:
                self.problem_solution = False
                print('Error: podczas rozwiązywania problemu P-median wystąpił błąd solvera Pulp - sprawdź poprawność wprowadzonych danych')

        elif self.algorytm == 'P-center':
            pcen_from_cm = spopt.locate.PCenter.from_cost_matrix(
                cost_matrix=self.cosMatrix,
                name=f"p-median-BikeStations-{self.obszar.name}",
                p_facilities=self.solutionNumber,
                predefined_facilities_arr=predef
            )
            try:
                self.problem_solution = pcen_from_cm.solve(solver)
            except:
                self.problem_solution = False
                print('Error: Podczas rozwiązywania problemu P-median wystąpił błąd solvera Pulp - sprawdź poprawność wprowadzonych danych')
        return self.problem_solution
    def createCostsMatrix(self):
        ntw = spaghetti.Network(in_data=self.net.NetDataFrame)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # ignore deprecation warning - GH pysal/libpysal#468
            ntw.snapobservations(self.pois, "clients", attribute=True)
        clients_snapped = spaghetti.element_as_gdf(ntw, pp_name="clients", snapped=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # ignore deprecation warning - GH pysal/libpysal#468
            ntw.snapobservations(self.propozycje, "facilities", attribute=True)
        facilities_snapped = spaghetti.element_as_gdf(ntw, pp_name="facilities", snapped=True)
        cost_matrix = ntw.allneighbordistances(
            sourcepattern=ntw.pointpatterns["clients"],
            destpattern=ntw.pointpatterns["facilities"],
        )
        print(cost_matrix.shape)
        return cost_matrix

    def getResults(self, resultsPath):
        print(datetime.now(), 'zapisuje rozwiazanie')
        self.solutions = gpd.GeoDataFrame(columns=['geometry', 'weights_sum', "how_many_points",'exists'])
        self.solutions.set_geometry('geometry')
        tmp_connections = []
        tmp_served = []
        if self.problem_solution == False:
            print('Wystąpił błąd podczas próby rozwiązania problemu. Nie można zwrócic wyniku')
            return 1
        for i, dv in enumerate(self.problem_solution.fac_vars):
            if dv.varValue:  # okresla czy jest w rozwiazaniu - 1 jest
                geom, predef = self.propozycje.iloc[i]
                ## jakie punkty są obsługiwane obslugiwane
                points = self.problem_solution.fac2cli[i]
                nr_points = len(points)
                weights_sum = 0
                if nr_points != 0:
                #konkretne punkty które są połączone z tym własnie rozwiązaniem
                    weights = self.weights.iloc[points]
                    for weight_idx, obsugiwany in enumerate(self.pois.iloc[points].geometry):
                        tmp_connections.append(shapely.LineString([geom, obsugiwany]))
                        tmp_served.append(obsugiwany)
                        weights_sum += weights.iloc[weight_idx]

                self.solutions.loc[i] = [geom,  weights_sum, nr_points, predef]
        solutions_to_save = gpd.GeoDataFrame(self.solutions[["geometry", "how_many_points",'weights_sum','exists']])
        solutions_to_save.set_geometry('geometry')
        print('Percentage covered: ', self.problem_solution.perc_cov)
        solutions_to_save = solutions_to_save.loc[solutions_to_save["exists"].isin(
            [False])]
        solutions_to_save.to_file(f"{resultsPath}/{self.algorytm}{self.obszar.name}s{self.skok}r{self.rad}_solution.gml", crs="EPSG:2180")
        connections_to_save = gpd.GeoSeries(tmp_connections)
        connections_to_save.to_file(f"{resultsPath}/{self.algorytm}{self.obszar.name}s{self.skok}r{self.rad}_solution_connections.gml", crs="EPSG:2180")
        served_to_save = gpd.GeoSeries(tmp_served)
        served_to_save.to_file(f"{resultsPath}/{self.algorytm}{self.obszar.name}s{self.skok}r{self.rad}_served.gml", crs="EPSG:2180")
        self.existingStations.to_file(f"{resultsPath}/{self.algorytm}{self.obszar.name}s{self.skok}r{self.rad}_existing_stations.gml", crs="EPSG:2180")
        return self.solutions
