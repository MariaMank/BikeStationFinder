import overpass
import geopandas as gpd
import pandas as pd
import shapely
from bs4 import BeautifulSoup
import numpy as np
import requests
import pulp
import warnings
from datetime import datetime
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # ignore deprecation warning - GH pysal/spaghetti#649
    import spaghetti
    import spopt
    from spopt.locate import MCLP, simulated_geo_points

from contextlib import suppress
from geopy import Nominatim


class Obszar:
    def __init__(self, name, admin_level):

        geolocalizador = Nominatim(user_agent="Obszar")
        localizacion = geolocalizador.geocode(name)
        if localizacion is None:
            raise ValueError("Nie moglam znalezc obszaru")
        self.area_id = 3600000000 + localizacion.raw['osm_id']
        #self.dataInFrame = self.import_data(name, admin_level)

    def import_data(self, name, admin_level):
        api = overpass.API()
        result = api.get(f"""
                area({self.area_id})->.b;
                rel['boundary'='administrative']['admin_level'='{admin_level}'](area.b);
            """, verbosity='geom')
        print(result)
        dane1 = gpd.GeoDataFrame.from_features(result,geometry='geometry', crs="EPSG:4326")  # 3857
        # print(footways1["geometry"][0], footways1["geometry"][0].length)
        dane2 = dane1.to_crs(epsg=2180)
        print(dane2)
        # print(footways1["geometry"][0], footways1["geometry"][0].length)
        # footways2.to_file(f"./data/{key}.geojson", crs="EPSG:2180")
        return dane2

    def saving(self, name):
        self.dataInFrame.to_file(f"./noweDane/{name}.geojson", crs="EPSG:2180")
        print(f"zapisano GeodataFrame z drogami Residential do ./noweDane/{name}.geojson ")



mójObszar = Obszar('Bielany', 7)

class DataI():
    def __init__(self, query:str, obszar_id:int):
        self.dataInFrame = self.import_data(query, obszar_id)


    def import_data(self, query, obszar_id):
        api = overpass.API()
        try:
            #area({obszar_id})->.b;
            result = api.get(f"""
                    area['name'='Mińsk Mazowiecki']['admin_level'=8]->.b;
                    (
                    {query};
                );
                
                """, verbosity='geom',build=True)
            print(query)
        except overpass.errors.UnknownOverpassError:
            result = {'features':[]} #nie moze tak zostac
            print('ehh')
            #zaimplementować żeby parsowało i zwracało bez wadliwych obiektow - waye czescia relacji budynkow
        if (len(result['features']) == 0):
            dane1 = gpd.GeoDataFrame()
        else:
            dane1 = gpd.GeoDataFrame.from_features(result, crs="EPSG:4326")  # 3857
        # print(footways1["geometry"][0], footways1["geometry"][0].length)
            dane1 = dane1.to_crs(epsg=2180)
        # print(footways1["geometry"][0], footways1["geometry"][0].length)
        # footways2.to_file(f"./data/{key}.geojson", crs="EPSG:2180")
        return dane1

    def saving(self, name):
        self.dataInFrame.to_file(f"./noweDane/{name}.geojson", crs="EPSG:2180")
        print(f"zapisano GeodataFrame z drogami Residential do ./noweDane/{name}.geojson ")


                                                                                                                                                                                                                                                    #626505.9068687708 481238.6833222164 631139.0765120629 484272.46604443155

class Proposition:
    def __init__(self, skok, bariery, siec):
        self.siec = siec
        [minx1, miny1, maxx1, maxy1] = siec.total_bounds
        x_diff = maxx1 - minx1
        y_diff = maxy1 - miny1
        minx1 = minx1 + x_diff/3
        maxx1 = maxx1 - x_diff/3
        miny1 = miny1 + y_diff/3
        maxy1 = maxy1 - y_diff/3
        Ursus = self.get_UrsusPolygon(minx1, miny1, maxx1, maxy1) #geojson zbiór dzielnic, ale ma tylko jedna pozycje, bo bylo zapytanie o nazwe
        #Ursus = Ursus.to_crs(epsg=2180)
        [minx, miny, maxx, maxy] = Ursus.total_bounds
        minx = int(minx) + 1
        miny = int(miny) + 1
        maxx = int(maxx)
        maxy = int(maxy)
        self.siatka = self.ceratePoints(minx, miny, maxx, maxy, skok, Ursus[0]) #zeby sie nie tworzyly w kwadracie ale poza obszarem ursusa
        self.bariery = bariery.all_barriers
        self.dataInFrame = gpd.GeoDataFrame(self.filtration(), crs="EPSG:2180")
        self.dataInFrame.to_file(f"./noweDane/propozycjeP.geojson", crs="EPSG:2180")

    #przepisac na prg!!! - obwody nie cos innego
    def get_UrsusPolygon(self, minx, miny, maxx, maxy):
        #ursus_wfs = "http://mapa.um.warszawa.pl/WebServices/GraniceDzielnic/wgs84/findByName/Ursus"

        ursus_wfs = f"https://mapy.geoportal.gov.pl/wss/service/PZGIK/PRG/WFS/AdministrativeBoundaries/wfs?request=getFeature&version=2.0.0&service=WFS&typename=ms:A05_Granice_jednostek_ewidencyjnych&bbox={miny},{minx},{maxy},{maxx}"
        response = requests.get(ursus_wfs)
        #sparsowac to recznie :((
        geom, name = self.parsePRGgml(response)
        geom = shapely.Polygon(geom)
        dane1 = gpd.GeoSeries(geom, crs="EPSG:2180")
        dane1.to_file(f"./noweDane/{name}.geojson", crs="EPSG:2180")
        #dane1 = gpd.read_file(response.text, crs="EPSG:2180")]
        return dane1

    def parsePRGgml(self, response):
        gml = BeautifulSoup(response.text, features='xml')
        wszystkieGranice = gml.find_all('ms:A05_Granice_jednostek_ewidencyjnych')
        if len(wszystkieGranice) == 1:
            nazwa = wszystkieGranice[0].contents[11].contents[0]
            geom = wszystkieGranice[0].find_all('gml:posList')[0].contents[0].split()
            newGeom = []
            temp=[]
            for id, i in enumerate(geom):
                if id%2 == 1:
                    temp.append(i)
                    temp.append(y)
                    newGeom.append(shapely.Point(temp))
                    temp=[]
                else:
                    y=i
        return newGeom, nazwa

    def ceratePoints(self, minx, miny, maxx, maxy, skok, Ursus):
        tmp = []
        startx = minx
        while startx<maxx:
            starty = miny
            while starty< maxy:
                if shapely.within(shapely.Point([startx, starty]), Ursus):
                    newPoint = shapely.Point([startx, starty])
                    tmp.append(newPoint)
                starty += skok
            startx+=skok
        return gpd.GeoSeries(tmp).set_crs(crs="EPSG:2180")

    def filtration(self):
        inside = self.siatka.clip(self.bariery) #wybor punktow tylko w budynkach i drogach
        new_inside = gpd.GeoDataFrame({'geometry': inside}, crs="EPSG:2180")
        outside = gpd.GeoDataFrame({'geometry': self.siatka}, crs="EPSG:2180") #wszystkie punkty
        return outside.overlay(new_inside, "difference") #roznica tych dwoch zbiorow - czyli te punkty ktore nie sa ani w drogach ani w budynach


class Net:
    def __init__(self):
        residential_query = 'way["highway"="residential"]["sidewalk"!="yes"]["sidewalk"!="both"]["sidewalk"!="separate"]["sidewalk"!="left"]["sidewalk"!="right"]["sidewalk:left" != "separate"]["sidewalk:right" != "separate"](area.b)'
        residential = DataI(residential_query, mójObszar.area_id)
        pedestrians = DataI(self.__createHghwayQuery("pedestrian"), mójObszar.area_id)
        steps = DataI(self.__createHghwayQuery("steps"), mójObszar.area_id)
        serviceRoads = DataI(self.__createHghwayQuery("service"), mójObszar.area_id)
        livingStreet = DataI(self.__createHghwayQuery("living_street"), mójObszar.area_id)
        footways = DataI(self.__createHghwayQuery("footway"), mójObszar.area_id)
        paths_query = 'way["highway"="path"]["foot"="designated"](area.b);' \
                      'way["informal"="yes"](area.b)'
        paths = DataI(paths_query, mójObszar.area_id)
        self.NetDataFrame = pd.concat([residential.dataInFrame.geometry,pedestrians.dataInFrame.geometry, serviceRoads.dataInFrame.geometry, livingStreet.dataInFrame.geometry,
                                      footways.dataInFrame.geometry, paths.dataInFrame.geometry, steps.dataInFrame.geometry]).to_crs(crs="EPSG:2180")
        self.NetDataFrame = self.NetDataFrame.drop_duplicates()
        #print(self.NetDataFrame.iloc[:10])
        #self.NetDataFrame.to_file(f"./noweDane/NetP.geojson", crs="EPSG:2180")

    def __createHghwayQuery(self, name):
        return f'way["highway"="{name}"](area.b)'

    def createConnectionToTheStations(self, propositions, existing, budynki):
        newLines = gpd.GeoSeries()
        # print(idx_odleglosci)
        for i in propositions.geometry:
            odleglosci = self.NetDataFrame.geometry.distance(i)
            minDist = np.min(odleglosci)
            shortesLine = shapely.shortest_line(i, self.NetDataFrame.geometry.iloc[np.where(odleglosci== minDist)[0]])
            newLines = pd.concat([newLines.geometry, shortesLine.geometry])
        for i in existing.geometry:
            odleglosci = self.NetDataFrame.geometry.distance(i)
            minDist = np.min(odleglosci)
            shortesLine = shapely.shortest_line(i, self.NetDataFrame.geometry.iloc[np.where(odleglosci == minDist)[0]])
            newLines = pd.concat([newLines.geometry, shortesLine.geometry])
        for i in budynki.geometry:
            odleglosci = self.NetDataFrame.geometry.distance(i)
            minDist = np.min(odleglosci)
            shortesLine = shapely.shortest_line(i, self.NetDataFrame.geometry.iloc[np.where(odleglosci == minDist)[0]])
            newLines = pd.concat([newLines.geometry, shortesLine.geometry])
        #newLines.to_file(f"./noweDane/newLines.geojson", crs="EPSG:2180")
        self.NetDataFrame = gpd.GeoDataFrame(geometry=pd.concat([self.NetDataFrame.geometry, newLines.geometry]))
        self.NetDataFrame.set_crs(crs="EPSG:2180")
        self.NetDataFrame.to_file(f"./noweDane/NetP.geojson", crs="EPSG:2180")
        return self.NetDataFrame


class Bariers:
    def __init__(self):
        buildings = DataI('rel["building"](area.b)', mójObszar.area_id)
        buildingsP = self.buildingsToPolygons(buildings.dataInFrame.geometry)
        railways = DataI('way["railway" = "rail"](area.b)', mójObszar.area_id)
        ways1 = DataI('way["highway" = "motorway"](area.b); way["highway"= "primary"](area.b); way["highway" = "secondary"](area.b); way["highway"="tertiary"](area.b)', mójObszar.area_id)
        ways2 = DataI('way["highway"="service"](area.b);way ["highway"="living_street"](area.b); way["highway"="residential"](area.b);way["highway"="trunk"](area.b)', mójObszar.area_id)
        way1B = self.waysToPolygons(pd.concat([ways1.dataInFrame.geometry, railways.dataInFrame.geometry]).to_crs(crs="EPSG:2180"), 4)
        way2B = self.waysToPolygons(ways2.dataInFrame.geometry, 2)
        self.all_barriers = pd.concat([way1B, way2B, buildingsP]).to_crs(crs="EPSG:2180")
        #self.all_barriers.to_file(f"./noweDane/bariers.geojson", crs="EPSG:2180")

    def waysToPolygons(self, ways, buff):
        return ways.buffer(buff)

    def buildingsToPolygons(self, buildings_geom):
        tmp = []
        for index, poi in enumerate(buildings_geom):
            try:
                new_point = shapely.Polygon(buildings_geom.loc[index].coords)
                tmp.append( new_point)
            except:
                tmp.append(buildings_geom.loc[index])

        return gpd.GeoSeries(tmp).set_crs(crs="EPSG:2180")



###ZOSTALO


#POBRAC BUDYNKI - CENTROIDY - ZMIERZYC POWIERZCHNIE, NA TEJ PODSTAWIE POWAGOWAĆ - done

class Building:
    def __init__(self):
        self.buildings = DataI('nwr["building"](area.b)', mójObszar.area_id)
        self.new_buildings = self.buildings.dataInFrame[["building", "shop", "amenity", "level", "surface", "house", "building:levels", "area", "geometry"]]
        self.toCentroids()
        self.living_houses = self.new_buildings.loc[self.new_buildings["building"].isin([ 'detached', 'apartments','presbytery', 'house', 'residential' ,'yes','semidetached_house', 'centroid'])]
        self.find_area()
        #print("suma: ",self.living_houses["people"].describe())
        #self.living_houses.set_geometry("centroid",crs="EPSG:2180" )
        #nowy = gpd.GeoDataFrame({'geometry':self.living_houses["centroid"], 'total_area':self.living_houses["total_area"] }, crs="EPSG:2180")
        #nowy.to_file(f"./noweDane/living_houses.geojson", crs="EPSG:2180")
        #self.living_houses["centroid"].to_file(f"./noweDane/living_houses.geojson", crs="EPSG:2180")
        self.living_houses = gpd.GeoDataFrame({'geometry':self.living_houses["centroid"], 'weights':self.living_houses["weights"] }, crs="EPSG:2180")

    def toCentroids(self):
        centroidy = self.new_buildings["geometry"].centroid
        self.new_buildings = self.new_buildings.assign(centroid= centroidy)

    def find_area(self):
        tmp = []
        wagi_tmp = []
        for index, poi in self.living_houses.iterrows():
            try:
                #new_area = shapely.Polygon(self.living_houses.iloc[index]["geometry"].coords).area
                new_area = shapely.Polygon(poi["geometry"].coords).area*int(poi["building:levels"])
                tmp.append(new_area)
            except ValueError:
                new_area=shapely.Polygon(poi["geometry"].coords).area
                tmp.append(new_area)
            except NotImplementedError:
                new_area = poi["geometry"].area*int(poi["building:levels"])
                tmp.append(new_area)
            #except IndexError:
            #    new_area = shapely.Polygon(poi["geometry"].coords)
            if new_area <= 200:
                wagi_tmp.append(2.5)
            elif new_area <=400:
                wagi_tmp.append(4)
            else:
                waga = new_area*0.85/30
                wagi_tmp.append(waga)
        #self.living_houses= self.living_houses.assign(total_area=tmp, people=wagi_tmp)
        self.living_houses = self.living_houses.assign(weights=wagi_tmp)


        #return gpd.GeoSeries(tmp).set_crs(crs="EPSG:2180")
    def give_weights(self):
        return self.living_houses["weights"]


class BikeStationSolver:
    def __init__(self, service_radius, solutionsNumber,frequencyPropPoints,area=None, existing=None, net=None, proposiotions=None, weights=None, pois=None):
        self.rad = service_radius
        self.skok = frequencyPropPoints
        print(datetime.now(), "start")
        print(datetime.now(),"tworze siec")
        if (net == None):
            self.net = Net()  # siec piesza
            #self.net = self.net.NetDataFrame
        else:
            self.net = net
        print(datetime.now(),'pobieram dane o istniejacych stacjach')
        if (existing == None):
            self.existingStations = DataI(
                "node[amenity = bicycle_rental](area.b)", mójObszar.area_id)  # istniejace stacje rowerow - kokurencja
            if (self.existingStations.dataInFrame.size != 0):
                self.existingStations = gpd.GeoDataFrame(self.existingStations.dataInFrame.geometry)
                self.existingStations = self.existingStations.assign(predefined_loc=np.ones(len(self.existingStations.geometry),  dtype=bool))
        else:
            self.existingStations = existing
        # existing_bikeSharing.saving("existingP")
        if (self.existingStations.dataInFrame.size != 0):
            self.solutionNumber = solutionsNumber + len(self.existingStations.geometry)
        print(datetime.now(),'tworze punkty na potencjalne stacje')
        self.bariery = Bariers()  # bariery - tam gdzie nie moze byc stacji - budynki i drogi, tory
        if (proposiotions != None):
            self.propozycje = proposiotions
        else:
            self.propozycje = Proposition(self.skok, self.bariery,self.net.NetDataFrame)  # kwadratowa siatka propozycji (punkty) po selekcji tego co juz nie moze byc - jakieś 1.5 min
            self.propozycje = self.propozycje.dataInFrame
            self.propozycje = self.propozycje.assign(predefined_loc=np.zeros(len(self.propozycje.geometry),  dtype=bool))
            self.propozycje = pd.concat([self.propozycje, self.existingStations])

        print(datetime.now(),'pobieram i waguje budynki mieszkaniowe')
        if (pois == None):
            self.budynki = Building()
            if (weights == None):
                self.weights = self.budynki.give_weights()
            else:
                self.weights = weights
            self.budynki = self.budynki.living_houses

        else:
            self.budynki = pois
        print(datetime.now(),'tworze polaczenia propozycji, pois i budynków mieszkaniowych z reszta sieci')
        self.net.createConnectionToTheStations(self.propozycje, self.existingStations, self.budynki)

    def solve(self):
        print(datetime.now(),'szukam rozwiazania')
        # set the solver
        solver = pulp.COIN_CMD(msg=False, warmStart=True, keepFiles=True, path=r'D:\OneDrive - Politechnika Warszawska\inżynierka\inz\env\Lib\site-packages\pulp\solverdir\cbc\win\64\cbc.exe')
        print(datetime.now(),'tworze macierz odleglosci po sieci')
        self.cosMatrix = self.createCostsMatrix()
        predef = self.propozycje['predefined_loc'].array
        mclp_from_cm = MCLP.from_cost_matrix(
            self.cosMatrix,
            self.weights,
            self.rad,
            p_facilities=self.solutionNumber,
            name="mclp-BikeStations-Legionowo",
            predefined_facilities_arr= predef
        )
        print(datetime.now(),'zaczynam rozwiazanie')
        self.mclp_from_cm = mclp_from_cm.solve(solver)
        return self.mclp_from_cm.fac_vars

    def createCostsMatrix(self):
        ntw = spaghetti.Network(in_data=self.net.NetDataFrame)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # ignore deprecation warning - GH pysal/libpysal#468
            ntw.snapobservations(self.budynki, "clients", attribute=True)
        clients_snapped = spaghetti.element_as_gdf(ntw, pp_name="clients", snapped=True)
        #clients_snapped.drop(columns=["id", "comp_label"], inplace=True)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # ignore deprecation warning - GH pysal/libpysal#468
            ntw.snapobservations(self.propozycje, "facilities", attribute=True)
        facilities_snapped = spaghetti.element_as_gdf(ntw, pp_name="facilities", snapped=True)
        #facilities_snapped.drop(columns=["id", "comp_label"], inplace=True)
        cost_matrix = ntw.allneighbordistances(
            sourcepattern=ntw.pointpatterns["clients"],
            destpattern=ntw.pointpatterns["facilities"],
        )
        print(cost_matrix, cost_matrix.shape)
        return cost_matrix
    def getResults(self):
        print(datetime.now(),'zapisuje rozwiazanie')
        self.solutions = gpd.GeoDataFrame(columns=['geometry', 'pointsServiced', "how_many_points", 'exists'])
        self.solutions.set_geometry('geometry')
        for i, dv in enumerate(self.mclp_from_cm.fac_vars):
            if dv.varValue:  # okresla czy jest w rozwiazaniu - 1 jest
                geom, predef = self.propozycje.iloc[i]
                ## what points are obslugiwane
                points = self.mclp_from_cm.fac2cli[i]
                # jakby chciec konkretne punkty
                # self.budynki.iloc[self.mclp_from_cm.fac2cli[i]]
                nr_points = len(points)
                self.solutions.loc[i] = [geom, points, nr_points, predef]
        to_save = gpd.GeoDataFrame(self.solutions[["geometry", "how_many_points", 'exists']])
        to_save.set_geometry('geometry')
        print('Percentage covered: ', self.mclp_from_cm.perc_cov)
        to_save.to_file(f"./noweDane/Solution.geojson", crs="EPSG:2180")
        return self.solutions

        # facility-(client) symbology and legend entries


solver = BikeStationSolver(500, 3, 200)
#napisac to dla istniejacych plikow
solver.solve()
solver.getResults()
print(solver.solutions)
print('ockm')

"""
#
# quantity demand points
CLIENT_COUNT = 100

# quantity supply points
FACILITY_COUNT = 10

# maximum service radius (in distance units)
SERVICE_RADIUS = 4

# number of candidate facilities in optimal solution
P_FACILITIES = 4

# random seeds for reproducibility
CLIENT_SEED = 5
FACILITY_SEED = 6
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # ignore deprecation warning - GH pysal/libpysal#468
    lattice = spaghetti.regular_lattice((0, 0, 10, 10), 9, exterior=True)
ntw = spaghetti.Network(in_data=lattice)
streets = spaghetti.element_as_gdf(ntw, arcs=True)
streets_buffered = gpd.GeoDataFrame(
    gpd.GeoSeries(streets["geometry"].buffer(0.5).unary_union),
    crs=streets.crs,
    columns=["geometry"],
)
client_points = simulated_geo_points(
    streets_buffered, needed=CLIENT_COUNT, seed=CLIENT_SEED
)
facility_points = simulated_geo_points(
    streets_buffered, needed=FACILITY_COUNT, seed=FACILITY_SEED
)
np.random.seed(0)
ai = np.random.randint(1, 12, CLIENT_COUNT)
client_points["weights"] = ai

#solver = BikeStationSolver(500, 3, 200, net=streets, proposiotions=facility_points, pois=client_points, weights=ai )
def plot_results(model, p, facs, clis=None, ax=None):
 """#Visualize optimal solution sets and context."""
"""

    # extract facility-client relationships for plotting (except for p-dispersion)
    fac_sites = {}
    for i, dv in enumerate(model.fac_vars):
        if dv.varValue:     #okresla czy jest w rozwiazaniu - 1 jest
            dv, predef = facs.loc[i, ["dv", "predefined_loc"]]
            fac_sites[dv] = [i, predef]
            #what points are obslugiwane
            geom = clis.iloc[model.fac2cli[i]]["geometry"]
            cli_points[dv] = geom

    # facility-(client) symbology and legend entries
    zorder = 4
    for fname, (fac, predef) in fac_sites.items(): #fname = dw, fac = index of fac
        cset = dv_colors[fname]
        if plot_clis:
            # clients
            geoms = cli_points[fname]
            gdf = geopandas.GeoDataFrame(geoms)
            gdf.plot(ax=ax, zorder=zorder, ec="k", fc=cset, markersize=100 * markersize)
            _label = f"Demand sites covered by {fname}"
            _mkws = dict(markerfacecolor=cset, markeredgecolor="k", ms=markersize + 7)
            legend_elements.append(
                mlines.Line2D([], [], marker="o", lw=0, label=_label, **_mkws)
            )
        # facilities
        ec = "k"
        lw = 2
        predef_label = "predefined"
        if model.name.endswith(predef_label) and predef:
            ec = "r"
            lw = 3
            fname += f" ({predef_label})"
        facs.iloc[[fac]]



#print(budynki.buildings.dataInFrame.columns)



#POBRAĆ POI!!!!
#NADAĆ WAAGI


#PUŚCIC CAŁY ALGORYTM


#poprawic - geocoding, pobieranie z PRG,

plot_results(mclp_from_cm, P_FACILITIES, facility_points, clis=client_points)
def plot_results(model, p, facs, clis=None, ax=None):
    #Visualize optimal solution sets and context.
    if not ax:
        multi_plot = False
        fig, ax = plt.subplots(figsize=(6, 6))
        markersize, markersize_factor = 4, 4
    else:
        ax.axis("off")
        multi_plot = True
        markersize, markersize_factor = 2, 2
    ax.set_title(model.name, fontsize=15)

    # extract facility-client relationships for plotting (except for p-dispersion)
    plot_clis = isinstance(clis, geopandas.GeoDataFrame)
    if plot_clis:
        cli_points = {}
    fac_sites = {}
    for i, dv in enumerate(model.fac_vars):
        if dv.varValue:
            dv, predef = facs.loc[i, ["dv", "predefined_loc"]]
            fac_sites[dv] = [i, predef]
            if plot_clis:
                geom = clis.iloc[model.fac2cli[i]]["geometry"]
                cli_points[dv] = geom

    # study area and legend entries initialization
    streets.plot(ax=ax, alpha=1, color="black", zorder=1)
    legend_elements = [mlines.Line2D([], [], color="black", label="streets")]

    if plot_clis:
        # any clients that not asscociated with a facility
        if model.name.startswith("mclp"):
            c = "k"
            if model.n_cli_uncov:
                idx = [i for i, v in enumerate(model.cli2fac) if len(v) == 0]
                pnt_kws = dict(ax=ax, fc=c, ec=c, marker="s", markersize=7, zorder=2)
                clis.iloc[idx].plot(**pnt_kws)
            _label = f"Demand sites not covered ($n$={model.n_cli_uncov})"
            _mkws = dict(marker="s", markerfacecolor=c, markeredgecolor=c, linewidth=0)
            legend_elements.append(mlines.Line2D([], [], ms=3, label=_label, **_mkws))

    # all candidate facilities
    facs.plot(ax=ax, fc="brown", marker="*", markersize=80, zorder=8)
    _label = f"Facility sites ($n$={len(model.fac_vars)})"
    _mkws = dict(marker="*", markerfacecolor="brown", markeredgecolor="brown")
    legend_elements.append(mlines.Line2D([], [], ms=7, lw=0, label=_label, **_mkws))

    # facility-(client) symbology and legend entries
    zorder = 4
    for fname, (fac, predef) in fac_sites.items():
        cset = dv_colors[fname]
        if plot_clis:
            # clients
            geoms = cli_points[fname]
            gdf = geopandas.GeoDataFrame(geoms)
            gdf.plot(ax=ax, zorder=zorder, ec="k", fc=cset, markersize=100 * markersize)
            _label = f"Demand sites covered by {fname}"
            _mkws = dict(markerfacecolor=cset, markeredgecolor="k", ms=markersize + 7)
            legend_elements.append(
                mlines.Line2D([], [], marker="o", lw=0, label=_label, **_mkws)
            )
        # facilities
        ec = "k"
        lw = 2
        predef_label = "predefined"
        if model.name.endswith(predef_label) and predef:
            ec = "r"
            lw = 3
            fname += f" ({predef_label})"
        facs.iloc[[fac]].plot(
            ax=ax, marker="*", markersize=1000, zorder=9, fc=cset, ec=ec, lw=lw
        )
        _mkws = dict(markerfacecolor=cset, markeredgecolor=ec, markeredgewidth=lw)
        legend_elements.append(
            mlines.Line2D([], [], marker="*", ms=20, lw=0, label=fname, **_mkws)
        )
        # increment zorder up and markersize down for stacked client symbology
        zorder += 1
        if plot_clis:
            markersize -= markersize_factor / p

    if not multi_plot:
        # legend
        kws = dict(loc="upper left", bbox_to_anchor=(1.05, 0.7))
        plt.legend(handles=legend_elements, **kws)


"""
#import processing
#processing.runalg('qgis:distancetonearesthub', points, hubs, field, geometry, unit, output)