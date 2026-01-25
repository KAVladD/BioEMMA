import pathlib
from urllib.request import urlretrieve

from bs4 import BeautifulSoup

class KeggMap():

    def __init__(self, 
                 metabolites: list = [], 
                 metabolites_positions: dict = {},
                 reactions: list = [], 
                 reaction_positions: dict = {}
                 ) -> None:

        self.reset()

        self.metabolites = metabolites
        self.metabolites_positions = metabolites_positions
        self.reactions = reactions
        self.reaction_positions = reaction_positions 

    def read_from_file(self, file):

        self.reset()

        kegg_xml = BeautifulSoup(file, features="xml")

        self._get_metabolites_from_xml(kegg_xml)
        self._get_reactions_from_xml(kegg_xml)

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

    def reset(self):

        self.metabolites = []
        self.metabolites_positions = {}
        self.reactions = []
        self.reaction_positions = {} 



if __name__ == "__main__":

    test_url = "https://rest.kegg.jp/get/rn00630/kgml"
    test_file = "glycoxylate.xml"

    m = KeggMap()

    # with open(test_file, "r") as f:
    #     m.read_from_file(f)

    m.read_from_url(test_url)



    #https://www.kegg.jp/kegg-bin/download?entry=rn00010&format=kgml
        