from maps import KeggMap
import numpy as np

from utils.id_dicts import make_metabolites_id_dict, make_reactions_id_dict


class EscherMapper:

    def __init__(self, cobra_model):

        self.cobra_model = cobra_model
        self._extend_metabolite_description()
        self._extend_reactions_description()

        self.markers_dist = 10
        self.factor = 4
        self.segments_counter = 0

        self.metabolite_label_shift = [10, 10]
        self.reaction_label_shift = [10, -20]

        self.visible_r_indxs = []
        self.visible_m_indxs = []
        self.visible_m_ids = []
        self.m_id_to_kegg_name = {}
        self.metabolite_id2node = {}
        self.reactions_idx2kegg = {}

        self.reactions = {}
        self.nodes = {}
        self.text_labels = {}
        self.canvas = {
            "x": 0,
            "y": 0,
            "width": 100,
            "height": 100,
        }

    def map_flux_on_kegg(self, kegg_map):

        escher_map = []

        description = self._generate_description()
        escher_map.append(description)

        visible_reactions = self._intersect_reactions(kegg_map.reactions)
        visible_metabolites = self._extract_metabolites(kegg_map.metabolites)

        # print(set(visible_metabolites))
        # print(set(visible_reactions))

        self._add_metbolites(kegg_map.metabolites_positions)
        self._add_reactions(kegg_map.reaction_positions)

        self._multiply_positions()
        self._tune_canvas()

        model = {}
        model["reactions"] = self.reactions
        model["nodes"] = self.nodes
        model["text_labels"] = self.text_labels
        model["canvas"] = self.canvas

        escher_map.append(model)

        return escher_map

    def _extend_metabolite_description(self):

        id_dict = make_metabolites_id_dict()

        for m in self.cobra_model.metabolites:

            bigg = m.id[:-2]
            seed = list(m.id.split("_"))[0]

            m.annotation["bigg.metabolite"] = bigg

            kegg = id_dict.get(bigg, None)

            if not kegg:
                kegg = id_dict.get(seed, None)

            #print(kegg, seed)
            if not kegg:
                continue


            m.annotation["kegg.compound"] = kegg

    def _extend_reactions_description(self):

        id_dict = make_reactions_id_dict()

        for r in self.cobra_model.reactions:

            bigg = r.id
            # print("kegg.reaction" in r.annotation.keys())

            try:
                kegg = id_dict[bigg]
                r.annotation["bigg.reaction"] = bigg
                r.annotation["kegg.reaction"] = kegg
            except:
                # print(bigg)
                pass
            # if m.id[-2:] == "_e":
            #     continue

            # kegg = id_dict[bigg]

            # m.annotation["bigg.metabolite"] = bigg
            # m.annotation["kegg.compound"] = kegg

    def _generate_description(self):

        descr = {}
        descr["map_name"] = str(self.cobra_model)
        descr["map_id"] = "p+p"  # idk how to create uid
        descr["map_description"] = ""
        descr["homepage"] = "https://escher.github.io"
        descr["schema"] = "https://escher.github.io/escher/jsonschema/1-0-0#"

        return descr

    def _extract_metabolites(self, metabolites):

        vis_ms = []

        for i in self.visible_r_indxs:

            r = self.cobra_model.reactions[i]
            r_ms = r.metabolites

            for m in r_ms:

                if "kegg.compound" in m.annotation.keys():
                    # print("kegg")
                    kegg_name = m.annotation["kegg.compound"]

                    if type(kegg_name) is list:
                        for name in kegg_name:
                            if name in metabolites:
                                vis_ms.append(name)
                                # print("added to list")
                    else:
                        if kegg_name in metabolites:
                            vis_ms.append(kegg_name)
                            # print("added to list")

        self._extract_vis_m_indxs(vis_ms)

        return vis_ms

    def _extract_vis_m_indxs(self, visible_metabolites):

        for i, m in enumerate(self.cobra_model.metabolites):
            if "kegg.compound" in m.annotation.keys():
                kegg_name = m.annotation["kegg.compound"]

                if type(kegg_name) is list:
                    for name in kegg_name:
                        if name in visible_metabolites:
                            self.visible_m_indxs.append(i)
                            self.visible_m_ids.append(m.id)
                            self.m_id_to_kegg_name[m.id] = name
                else:
                    if kegg_name in visible_metabolites:
                        self.visible_m_indxs.append(i)
                        self.visible_m_ids.append(m.id)
                        self.m_id_to_kegg_name[m.id] = kegg_name

    def _intersect_reactions(self, reactions):

        vis_rs = []

        for i, r in enumerate(self.cobra_model.reactions):
            if "kegg.reaction" in r.annotation.keys():

                kegg_names = r.annotation["kegg.reaction"]
                if not (type(kegg_names) is list):
                    kegg_names = [kegg_names]

                for kegg_name in kegg_names:
                    if kegg_name in reactions:
                        vis_rs.append(kegg_name)
                        self.visible_r_indxs.append(i)
                        self.reactions_idx2kegg[i] = kegg_name

        return vis_rs

    def _add_metbolites(self, m_positions):

        for i, m in enumerate(self.cobra_model.metabolites):

            # print(m)

            #m_dict = self._generate_metabolite_dict(m.annotation["bigg.metabolite"])
            m_dict = self._generate_metabolite_dict(m.id, m.name)
            if m.id[-1] == "e":
                continue
            # print(m.id)

            if i in self.visible_m_indxs:

                kegg_name = m.annotation["kegg.compound"]

                m_dict = self._update_metabolite_pos(
                    m_dict, m_positions[self.m_id_to_kegg_name[m.id]]
                )

                node_n = str(len(self.nodes))
                self.nodes[node_n] = m_dict
                self.metabolite_id2node[m.id] = node_n

    def _generate_metabolite_dict(self, bigg_id, name=None, primary=True):

        name = name if name else bigg_id
        m_dict = {
            "node_type": "metabolite",
            "x": None,
            "y": None,
            "bigg_id": bigg_id,
            "name": name,
            "label_x": None,
            "label_y": None,
            "node_is_primary": primary,
        }

        return m_dict

    def _update_metabolite_pos(self, metabolite_dict, pos):

        metabolite_dict["x"] = pos[0]
        metabolite_dict["y"] = pos[1]
        metabolite_dict["label_x"] = str(float(pos[0]) + self.metabolite_label_shift[0])
        metabolite_dict["label_y"] = str(float(pos[1]) + self.metabolite_label_shift[0])

        return metabolite_dict

    def _add_reactions(self, reaction_positions):

        for i in self.visible_r_indxs:

            r = self.cobra_model.reactions[i]
            # kegg_name = r.annotation["kegg.reaction"]
            kegg_name = self.reactions_idx2kegg[i]

            r_dict = self._genarate_reaction_dict(r)
            r_dict = self._update_mdict_metabolites(r_dict, r)

            r_dict["label_x"] = str(float(reaction_positions[kegg_name][0]) + self.reaction_label_shift[0])
            r_dict["label_y"] = str(float(reaction_positions[kegg_name][1]) + self.reaction_label_shift[1])

            general_reaction_node_n = self._add_reaction_node(
                reaction_positions[kegg_name]
            )
            reaction_nodes = self._compose_reaction_nodes(r, general_reaction_node_n)

            r_dict = self._add_edges_from_nodes(r_dict, reaction_nodes)
            self.reactions[str(len(self.reactions))] = r_dict

    def _genarate_reaction_dict(self, reaction):

        reaction_dict = {
            "name": reaction.id,
            "bigg_id": reaction.id,
            "reversibility": reaction.reversibility,
            "label_x": None,
            "label_y": None,
            "gene_reaction_rule": "",
            "genes": [],
            "metabolites": [],
            "segments": {},
        }

        return reaction_dict

    def _update_mdict_metabolites(self, r_dict, r):

        for m, coef in r.metabolites.items():
            # simple_metabolite_dict = {
            #     "bigg_id": m.annotation["bigg.metabolite"],
            #     "coefficient": coef,
            # }
            simple_metabolite_dict = {
                "bigg_id": m.id,
                "coefficient": coef,
            }
            r_dict["metabolites"].append(simple_metabolite_dict)

        return r_dict

    def _add_reaction_node(self, pos):

        r_node = self._generate_not_primary_node_dict("midmarker", pos)

        node_n = str(len(self.nodes))
        self.nodes[node_n] = r_node

        return node_n

    def _add_marker_node(self, pos):

        r_node = self._generate_not_primary_node_dict("multimarker", pos)

        node_n = str(len(self.nodes))
        self.nodes[node_n] = r_node

        return node_n

    def _generate_not_primary_node_dict(self, type, pos):

        node = {"node_type": type, "x": pos[0], "y": pos[1]}

        return node

    def _compose_reaction_nodes(self, reaction, reaction_node_n):

        nodes = [[], [], []]

        primary_nodes = [[], [reaction_node_n], []]
        not_primary_nodes = [[], [], []]

        for m, coef in reaction.metabolites.items():

            primary = True if m.id in self.visible_m_ids else False

            if coef > 0:
                if primary:
                    primary_nodes[2].append(self.metabolite_id2node[m.id])
                    nodes[2].append(self.metabolite_id2node[m.id])
                else:
                    #not_primary_nodes[2].append(m.annotation["bigg.metabolite"])
                    not_primary_nodes[2].append([m.id, m.name])
            else:
                if primary:
                    primary_nodes[0].append(self.metabolite_id2node[m.id])
                    nodes[0].append(self.metabolite_id2node[m.id])
                else:
                    # not_primary_nodes[0].append(m.annotation["bigg.metabolite"])
                    not_primary_nodes[0].append([m.id, m.name])

        mean_positions = self._calc_mean_position(primary_nodes)
        markers_vecs = self._calc_markers_vecs(mean_positions)

        marker1_pos = mean_positions[1] + markers_vecs[0] * self.markers_dist
        marker2_pos = mean_positions[1] + markers_vecs[1] * self.markers_dist

        marker1_node_n = self._add_marker_node(marker1_pos)
        marker2_node_n = self._add_marker_node(marker2_pos)

        nodes[1].append(marker1_node_n)
        nodes[1].append(reaction_node_n)
        nodes[1].append(marker2_node_n)

        for i in [0, 2]:
            for j, m in enumerate(not_primary_nodes[i]):

                m_dict = self._generate_metabolite_dict(*m, primary=False)

                n = 0 if i == 0 else 1
                node_pos = (
                    mean_positions[1]
                    + 2 * markers_vecs[n] * self.markers_dist
                    + (-((len(not_primary_nodes[i]) - 1) // 2) - 0.5 + j)
                    * 3
                    * self.markers_dist
                    * np.array([[0, -1], [1, 0]])
                    @ markers_vecs[0]
                )
                m_dict = self._update_metabolite_pos(m_dict, node_pos)

                node_n = str(len(self.nodes))
                self.nodes[node_n] = m_dict

                nodes[i].append(node_n)

        # TODO func to add all not primary nodes to jsons

        return nodes

    def _calc_mean_position(self, pr_nodes):

        mean_pos = []

        for nodes in pr_nodes:
            s = [0, 0]
            for node_n in nodes:
                s[0] += float(self.nodes[node_n]["x"])
                s[1] += float(self.nodes[node_n]["y"])

            denom = len(nodes) if nodes else 1

            s[0], s[1] = s[0] / denom, s[1] / denom
            mean_pos.append(s)

        return np.array(mean_pos)

    def _calc_markers_vecs(self, mean_poss):

        vec1 = mean_poss[0] - mean_poss[1]
        vec2 = mean_poss[2] - mean_poss[1]

        vec1 = vec1 / np.linalg.norm(vec1)
        vec2 = vec2 / np.linalg.norm(vec2)

        return [vec1, vec2]

    def _add_edges_from_nodes(self, r_dict, nodes):

        for node_n in nodes[0]:

            self.segments_counter += 1
            edge_dict = self._generate_edge_dict(node_n, nodes[1][0])
            r_dict["segments"][str(self.segments_counter)] = edge_dict

        self.segments_counter += 1
        edge_dict = self._generate_edge_dict(nodes[1][0], nodes[1][1])
        r_dict["segments"][str(self.segments_counter)] = edge_dict

        self.segments_counter += 1
        edge_dict = self._generate_edge_dict(nodes[1][1], nodes[1][2])
        r_dict["segments"][str(self.segments_counter)] = edge_dict

        for node_n in nodes[2]:

            self.segments_counter += 1
            edge_dict = self._generate_edge_dict(nodes[1][2], node_n)
            r_dict["segments"][str(self.segments_counter)] = edge_dict

        return r_dict

    def _generate_edge_dict(self, from_node, to_node):

        edge_dict = {
            "from_node_id": from_node,
            "to_node_id": to_node,
            "b1": None,
            "b2": None,
        }

        return edge_dict

    def _tune_canvas(self):

        x, y = 0, 0

        for i, compound_dict in self.nodes.items():

            if compound_dict["x"]:

                if x < float(compound_dict["x"]):
                    x = float(compound_dict["x"])
                if y < float(compound_dict["y"]):
                    y = float(compound_dict["y"])

        self.canvas["width"] = x * 1.1
        self.canvas["height"] = y * 1.1

    def _multiply_positions(self):

        for node in self.nodes.values():

            node["x"] = float(node["x"]) * self.factor
            node["y"] = float(node["y"]) * self.factor

            if "label_x" in node.keys():

                node["label_x"] = float(node["label_x"]) * self.factor
                node["label_y"] = float(node["label_y"]) * self.factor

        for reaction in self.reactions.values():

            reaction["label_x"] = float(reaction["label_x"]) * self.factor
            reaction["label_y"] = float(reaction["label_y"]) * self.factor
