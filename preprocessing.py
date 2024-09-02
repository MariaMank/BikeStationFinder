from baseClass import BaseClass
import requests
import geopandas as gpd
import shapely
from bs4 import BeautifulSoup

from geopy import Nominatim
import pandas as pd
import numpy as np



class Obszar():
    def __init__(self, name):
        geolocalizador = Nominatim(user_agent="Obszar")
        localizacion = geolocalizador.geocode(name)
        if localizacion is None:
            raise ValueError("Nie moglem znalezc obszaru")
        self.area_id = 3600000000 + localizacion.raw['osm_id']
        self.name = name

    def saving(self, name):
        self.dataInFrame.to_file(f"./dane/{name}.geojson", crs="EPSG:2180")

class Proposition:
    def __init__(self, skok, bariery, siec, ready=None):
        if (ready== None):
            self.siec = siec
            [minx1, miny1, maxx1, maxy1] = siec.total_bounds
            x_diff = maxx1 - minx1
            y_diff = maxy1 - miny1
            minx1 = minx1 + x_diff / 3
            maxx1 = maxx1 - x_diff / 3
            miny1 = miny1 + y_diff / 3
            maxy1 = maxy1 - y_diff / 3
            Ursus = self.get_UrsusPolygon(minx1, miny1, maxx1,
                                          maxy1)  # geojson zbiór dzielnic, ale ma tylko jedna pozycje, bo bylo zapytanie o nazwe
            # Ursus = Ursus.to_crs(epsg=2180)
            [minx, miny, maxx, maxy] = siec.total_bounds
            minx = int(minx) + 1
            miny = int(miny) + 1
            maxx = int(maxx)
            maxy = int(maxy)
            self.siatka = self.ceratePoints(minx, miny, maxx, maxy, skok,
                                            Ursus[0])  # zeby sie nie tworzyly w kwadracie ale poza obszarem ursusa
            self.bariery = bariery.all_barriers
            self.dataInFrame = gpd.GeoDataFrame(self.filtration(), crs="EPSG:2180")
        else:
            self.dataInFrame = gpd.read_file(ready)

    def get_UrsusPolygon(self, minx, miny, maxx, maxy):
        ursus_wfs = f"https://mapy.geoportal.gov.pl/wss/service/PZGIK/PRG/WFS/AdministrativeBoundaries/wfs?request=getFeature&version=2.0.0&service=WFS&typename=ms:A05_Granice_jednostek_ewidencyjnych&bbox={miny},{minx},{maxy},{maxx}"
        response = requests.get(ursus_wfs)
        geom, name = self.parsePRGgml(response)
        self.name=name
        geom = shapely.Polygon(geom)
        dane1 = gpd.GeoSeries(geom, crs="EPSG:2180")
        dane1.to_file(f"./noweDane/nowePodejscie/{name}.shp", crs="EPSG:2180")
        return dane1

    def parsePRGgml(self, response):
        gml = BeautifulSoup(response.text, features='xml')
        wszystkieGranice = gml.find_all('ms:A05_Granice_jednostek_ewidencyjnych')
        if len(wszystkieGranice) == 1:
            nazwa = wszystkieGranice[0].contents[11].contents[0]
            geom = wszystkieGranice[0].find_all('gml:posList')[0].contents[0].split()
            newGeom = []
            temp = []
            for id, i in enumerate(geom):
                if id % 2 == 1:
                    temp.append(i)
                    temp.append(y)
                    newGeom.append(shapely.Point(temp))
                    temp = []
                else:
                    y = i
        return newGeom, nazwa

    def ceratePoints(self, minx, miny, maxx, maxy, skok, obszar):
        tmp = []
        startx = minx
        while startx < maxx:
            starty = miny
            while starty < maxy:
                if shapely.within(shapely.Point([startx, starty]), obszar):
                    newPoint = shapely.Point([startx, starty])
                    tmp.append(newPoint)
                starty += skok
            startx += skok
        return gpd.GeoSeries(tmp).set_crs(crs="EPSG:2180")

    def filtration(self):
        inside = self.siatka.clip(self.bariery)  # wybor punktow tylko w budynkach i drogach
        new_inside = gpd.GeoDataFrame({'geometry': inside}, crs="EPSG:2180")
        outside = gpd.GeoDataFrame({'geometry': self.siatka}, crs="EPSG:2180")  # wszystkie punkty
        return outside.overlay(new_inside, "difference")  # roznica tych dwoch zbiorow - czyli te punkty ktore nie sa ani w drogach ani w budynach

class Net:
    def __init__(self, obszar: Obszar,ready= None):
        if ready== None:
            self.obszar = obszar
            residential_query = 'way["highway"="residential"]["sidewalk"!="yes"]["sidewalk"!="both"]["sidewalk"!="separate"]["sidewalk"!="left"]["sidewalk"!="right"]["sidewalk:left" != "separate"]["sidewalk:right" != "separate"](area.b)'
            residential = BaseClass(query=residential_query, obszar_id=obszar.area_id)
            pedestrians = BaseClass(query=self.__createHghwayQuery("pedestrian"), obszar_id=obszar.area_id)
            steps = BaseClass(query=self.__createHghwayQuery("steps"), obszar_id=obszar.area_id)
            serviceRoads = BaseClass(query=self.__createHghwayQuery("service"), obszar_id=obszar.area_id)
            livingStreet = BaseClass(query=self.__createHghwayQuery("living_street"), obszar_id=obszar.area_id)
            footways = BaseClass(query=self.__createHghwayQuery("footway"), obszar_id=obszar.area_id)
            paths_query = 'way["highway"="path"]["foot"="designated"](area.b);' \
                          'way["informal"="yes"](area.b)'
            paths = BaseClass(query=paths_query, obszar_id=obszar.area_id)
            self.NetDataFrame = pd.concat(
                [residential.dataInFrame.geometry, pedestrians.dataInFrame.geometry, serviceRoads.dataInFrame.geometry,
                 livingStreet.dataInFrame.geometry,
                 footways.dataInFrame.geometry, paths.dataInFrame.geometry, steps.dataInFrame.geometry]).to_crs(
                crs="EPSG:2180")
            self.NetDataFrame = self.NetDataFrame.drop_duplicates()
        else:
            self.NetDataFrame = gpd.read_file(ready)

    def __createHghwayQuery(self, name):
        return f'way["highway"="{name}"](area.b)'

    def createConnectionToTheStations(self, propositions, existing, budynki):
        newLines = gpd.GeoSeries()
        for i in propositions.geometry:
            odleglosci = self.NetDataFrame.geometry.distance(i)
            minDist = np.min(odleglosci)
            shortesLine = shapely.shortest_line(i, self.NetDataFrame.geometry.iloc[np.where(odleglosci == minDist)[0]])
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
        self.NetDataFrame = gpd.GeoDataFrame(geometry=pd.concat([self.NetDataFrame.geometry, newLines.geometry]))
        self.NetDataFrame.set_crs(crs="EPSG:2180")
        self.NetDataFrame.to_file(f"./noweDane/nowePodejscie/{self.obszar.name}_Siec.shp", crs="EPSG:2180")
        return self.NetDataFrame


class Bariers:
    def __init__(self, obszar: Obszar):
        buildings = BaseClass(query='nwr["building"](area.b)', obszar_id=obszar.area_id)
        buildingsP = self.buildingsToPolygons(buildings.dataInFrame)
        railways = BaseClass(query='way["railway" = "rail"](area.b)', obszar_id=obszar.area_id)
        ways1 = BaseClass(
            query='way["highway" = "motorway"](area.b); way["highway"= "primary"](area.b); way["highway" = "secondary"](area.b); way["highway"="tertiary"](area.b)',
            obszar_id=obszar.area_id)
        ways2 = BaseClass(
            query='way["highway"="service"](area.b);way ["highway"="living_street"](area.b); way["highway"="residential"](area.b);way["highway"="trunk"](area.b)',
            obszar_id=obszar.area_id)
        way1B = self.waysToPolygons(
            pd.concat([ways1.dataInFrame.geometry, railways.dataInFrame.geometry]).to_crs(crs="EPSG:2180"), 4)
        way2B = self.waysToPolygons(ways2.dataInFrame.geometry, 2)
        self.all_barriers = pd.concat([way1B, way2B, buildingsP]).to_crs(crs="EPSG:2180")

    def waysToPolygons(self, ways, buff):
        return ways.buffer(buff)

    def buildingsToPolygons(self, buildings):
        tmp = []
        for index, poi in buildings.iterrows():
            if poi["geometry"].geom_type != 'Point':
                try:
                    new_area = shapely.Polygon(poi["geometry"].coords)
                    tmp.append(new_area)
                except NotImplementedError:
                    new_area = poi["geometry"]
                    tmp.append(new_area)
        return gpd.GeoSeries(tmp).set_crs(crs="EPSG:2180")

class Building:
    def __init__(self, obszar: Obszar):
        self.buildings = BaseClass(query='nwr["building"](area.b)', obszar_id=obszar.area_id)
        self.new_buildings = self.buildings.dataInFrame[
            ["building", "building:levels", "geometry"]]
        self.toCentroids()
        self.living_houses = self.new_buildings.loc[self.new_buildings["building"].isin(
            ['detached', 'apartments', 'presbytery', 'house', 'residential', 'yes', 'semidetached_house', 'centroid'])]
        self.find_area()
        self.living_houses = gpd.GeoDataFrame(
            {'geometry': self.living_houses["centroid"], 'weights': self.living_houses["weights"]}, crs="EPSG:2180")

    def toCentroids(self):
        centroidy = self.new_buildings["geometry"].centroid
        self.new_buildings = self.new_buildings.assign(centroid=centroidy)

    def find_area(self):
        tmp = []
        wagi_tmp = []
        for index, poi in self.living_houses.iterrows():
            try:
                new_area = shapely.Polygon(poi["geometry"].coords).area * int(poi["building:levels"])
                tmp.append(new_area)
            except ValueError:
                new_area = shapely.Polygon(poi["geometry"].coords).area
                tmp.append(new_area)
            except NotImplementedError:
                try:
                    levels = int(poi["building:levels"])
                except:
                    levels=1
                new_area = poi["geometry"].area * levels
                tmp.append(new_area)
            if new_area <= 200:
                wagi_tmp.append(2.5)
            elif new_area <= 400:
                wagi_tmp.append(4)
            else:
                waga = new_area * 0.85 / 30
                if waga <= 400:
                    wagi_tmp.append(waga)
                else:
                    wagi_tmp.append(400)
        self.living_houses = self.living_houses.assign(weights=wagi_tmp)

    def give_weights(self):
        return self.living_houses["weights"]

class Pois:
    def __init__(self, obszar: Obszar, buildings: Building):
        self.living_buildings = buildings.living_houses
        self.other_buildings = buildings.new_buildings
        self.other_buildings = self.other_buildings.loc[self.other_buildings["building"].isin(
            ['retail', 'hotel', 'commercial', 'office', 'supermarket', 'religious', 'college',
             'hospital', 'kindergarten', 'museum', 'public', 'school', 'train_station', 'transportation', 'university',
             'sports_hall', 'sport_center'])]
        self.out_door = BaseClass(query='way["leisure" = "park"](area.b); way ["leisure"="pitch"](area.b)', obszar_id=obszar.area_id)
        self.out_door = self.out_door.dataInFrame
        centroidy = self.out_door["geometry"].centroid
        self.out_door = self.out_door.assign(centroid=centroidy)
        self.out_door = self.out_door.assign(weights=40)
        self.give_weights()
        self.join()

    def give_weights(self):
        idx_3 = self.other_buildings['building'].isin(['retail','supermarket', 'commercial', 'public'])
        idx_2 = self.other_buildings['building'].isin(['train_station', 'transportation'])
        idx_1 = self.other_buildings['building'].isin(['college', 'university'])
        idx_4 = self.other_buildings['building'].isin(['school', 'sport_center', 'sports_hall'])
        idx_5 = self.other_buildings['building'].isin(['office', 'religious', 'hospital'])
        idx_6 = self.other_buildings['building'].isin(['kindergarten', 'hotel'])
        weights_tmp = []
        for i in self.other_buildings.index:
            if idx_3[i]:
                weights_tmp.append(150)
                #weights_tmp.append(75)
            elif idx_2[i]:
                weights_tmp.append(400)
            elif idx_1[i]:
                weights_tmp.append(200)
            elif idx_4[i]:
                weights_tmp.append(80)
            elif idx_5[i]:
                weights_tmp.append(60)
            elif idx_6[i]:
                weights_tmp.append(40)
        self.other_buildings = self.other_buildings.assign(weights=weights_tmp)

        return self.other_buildings


    def join(self):
        centroidy_tmp = pd.concat([self.living_buildings.geometry, self.other_buildings['centroid'], self.out_door['centroid']])
        weights_tmp = pd.concat([self.living_buildings['weights'], self.other_buildings['weights'], self.out_door['weights']])

        self.pois_inDataFrame = gpd.GeoDataFrame(
            {'geometry': centroidy_tmp, 'weights': weights_tmp}, crs="EPSG:2180")

        return self.pois_inDataFrame

