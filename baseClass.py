import overpass
import geopandas as gpd
import pandas as pd
import os



class BaseClass():
    query: str
    def __init__(self, obszar_id: int, query=None,ready=None):
        self.obszar_id = obszar_id
        if query != None:
            self.query = query
        if ready == None:
            self.dataInFrame = self.import_data()
        else:
            self.dataInFrame = self.toDataFrame(ready)

    def import_data(self):
        api = overpass.API()
        try:
            result = api.get(f"""
                    area({self.obszar_id})->.b; 
                    (
                    {self.query};
                );

                """, verbosity='geom', build=True)
        except overpass.errors.UnknownOverpassError:
            result = {'features': []}  # nie moze tak zostac
            print('Błąd podczas pobierania danych:', query)
        if (len(result['features']) == 0):
            dane1 = gpd.GeoDataFrame({'geometry': [None]}, crs="EPSG:2180")
        else:
            dane1 = gpd.GeoDataFrame.from_features(result, crs="EPSG:4326")
            dane1 = dane1.to_crs(epsg=2180)
        return dane1

    def save(self, name:str, format:str):
        path = f"./dane/"
        isExist = os.path.exists(path)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path)
        self.dataInFrame.to_file(f"{path}{name}.{format}", crs="EPSG:2180")
        print(f"zapisano GeodataFrame do {path}{name}.{path}")

    def join(self, data:list, columns:list):
        dict = {}
        for i in columns:
            tmp = []
            for j in data:
                tmp.append(data[i])
            dict[i] = pd.concat(tmp)

        return gpd.GeoDataFrame(dict, crs="EPSG:2180")


    def toDataFrame(self, ready):
        pass