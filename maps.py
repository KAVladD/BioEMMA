import pathlib
from urllib.request import urlretrieve

from bs4 import BeautifulSoup
from bidict import bidict

from utils.metanetx_mapper import MetaNetXMapper

class KeggMap():

    def __init__(self, 
                 metabolites: list = [], 
                 metabolites_positions: dict = {},
                 reactions: list = [], 
                 reaction_positions: dict = {}
                 ) -> None:
        
        self.m_mapper = MetaNetXMapper("resources/metabolite_mapping.tsv", "first")
        self.r_mapper = MetaNetXMapper("resources/reaction_mapping.tsv", "first")

        self.reset()

        self.metabolites = metabolites
        self.metabolites_positions = metabolites_positions
        self.reactions = reactions
        self.reaction_positions = reaction_positions 

        self.reaction_substrates = {}  
        self.reaction_products = {}    
        self.reaction_types = {}   

    def read_from_file(self, file):

        self.reset()

        kegg_xml = BeautifulSoup(file, features="xml")

        self._get_metabolites_from_xml(kegg_xml)
        self._get_reactions_from_xml(kegg_xml)
        self._get_reaction_details_from_xml(kegg_xml)

        self.r_k2b, self.r_k2s = self._get_reactions_annotations("00010")

    def read_from_url(self, url):

        urlretrieve(url, "temp.xml")

        with open("temp.xml", "r") as f:
            self.read_from_file(f)

    def _get_metabolites_from_xml(self, xml):

        entrys = xml.find_all("entry", type="compound")
        for entry in entrys:
            entry_graphics = entry.graphics
            metabolite_name = entry_graphics["name"]
            metabolite_pos = (entry_graphics["x"],
                                entry_graphics["y"])

            if not metabolite_name in self.metabolites:
                self.metabolites.append(metabolite_name)
            self.metabolites_positions[metabolite_name] = metabolite_pos

    #TODO process multiple positions       
    def _get_reactions_from_xml(self, xml):
        entrys = xml.find_all("entry", type="reaction")
        for entry in entrys:
            entry_graphics = entry.graphics

            reactions = entry["reaction"]

            for reaction in reactions.split():

                try:

                    reaction_name = reaction.split(":")[1]
                    reaction_pos = (entry_graphics["x"],
                                    entry_graphics["y"])
                    
                    if reaction_name not in self.reactions:
                        self.reactions.append(reaction_name)
                    self.reaction_positions[reaction_name] = reaction_pos
                except:
                    pass

    def _get_reaction_details_from_xml(self, xml):

        metabolite_ids = set()
        for metab in self.metabolites:
            metabolite_ids.add(metab)
        
        reaction_elements = xml.find_all("reaction")
        
        for reaction_elem in reaction_elements:
            reaction_name = reaction_elem.get("name", "")

            if ":" in reaction_name:
                reaction_name = reaction_name.split(":")[1]
            
            if reaction_name not in self.reactions:
                continue
            
            reaction_type = reaction_elem.get("type", "unknown")
            
            substrates = []
            for substrate in reaction_elem.find_all("substrate"):
                cpd_name = substrate.get("name", "").replace("cpd:", "")
                if cpd_name:
                    substrates.append(cpd_name)
            
            products = []
            for product in reaction_elem.find_all("product"):
                cpd_name = product.get("name", "").replace("cpd:", "")
                if cpd_name:
                    products.append(cpd_name)
            
            main_substrates = [s for s in substrates if s in metabolite_ids]
            side_substrates = [s for s in substrates if s not in metabolite_ids]
            
            main_products = [p for p in products if p in metabolite_ids]
            side_products = [p for p in products if p not in metabolite_ids]
            
            self.reaction_substrates[reaction_name] = {
                'main': main_substrates,
                'side': side_substrates
            }
            
            self.reaction_products[reaction_name] = {
                'main': main_products,
                'side': side_products
            }
            
            self.reaction_types[reaction_name] = reaction_type

    def reset(self):

        self.metabolites = []
        self.metabolites_positions = {}
        self.reactions = []
        self.reaction_positions = {} 

        self.reaction_substrates = {}  
        self.reaction_products = {}    
        self.reaction_types = {} 

    def get_metabolites(self):

        metabolites = {
            me: {
                "ids": {
                        "KEGG": me, 
                        "BIGG": self.m_mapper.get(me).bigg if self.m_mapper.get(me) else None, 
                        "SEED": self.m_mapper.get(me).seed if self.m_mapper.get(me) else None
                       },
                "position": self.metabolites_positions[me],
            }
            for me in self.metabolites
        }  

        return metabolites
    
    def get_reactions(self):

        reactions = {
            r: {
                "ids": {
                        "KEGG": r, 
                        "BIGG": self.r_mapper.get(r).bigg if self.r_mapper.get(r) else None, 
                        "SEED": self.r_mapper.get(r).seed if self.r_mapper.get(r) else None
                       },
                "position": self.reaction_positions[r],
                "substrates": self.reaction_substrates.get(r, {"main": [], "side": []}),
                "products": self.reaction_products.get(r, {"main": [], "side": []}),
                "reversibility": self.reaction_types.get(r, None),
            }
            for r in self.reactions
        }

        print(self.r_k2b, self.r_k2s)

        return reactions
    
    # annotations

    def _get_reactions_annotations(self, kegg_map_id):

        with open(f"resources/BIGG/map{kegg_map_id}.tsv") as f:

            data = f.readlines()
            kegg2bigg = {}

            for line in data[1:]:

                if not line:
                    continue

                ids = line.split()
                kegg = ids[1].split(";")
                bigg = ids[0]

                kegg2bigg |= {k: bigg for k in kegg}

        with open(f"resources/SEED/map{kegg_map_id}.tsv") as f:

            data = f.readlines()
            kegg2seed = {}

            for line in data[1:]:

                if not line:
                    continue

                ids = line.split()
                if len(ids) == 2:
                    kegg = ids[1].split(";")
                else:
                    kegg = ids[2].split(";")
                seed = ids[0]

                kegg2seed |= {k: seed for k in kegg}

        return kegg2bigg, kegg2seed


if __name__ == "__main__":

    test_url = "https://rest.kegg.jp/get/rn00630/kgml"
    test_file = "glycoxylate.xml"

    m = KeggMap()

    # with open(test_file, "r") as f:
    #     m.read_from_file(f)

    m.read_from_url(test_url)



    #https://www.kegg.jp/kegg-bin/download?entry=rn00010&format=kgml
        