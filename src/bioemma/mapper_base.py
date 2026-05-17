import numpy as np

from bioemma.metanetx_mapper import MetaNetXMapper
from bioemma._resources import resource_path

class EscherMapper:

    def __init__(self, 
                 metabolites: dict, 
                 reactions: dict,
                 markers_dist: int = 10,
                 scaling_factor: float = 4,
                 metabolite_label_shift: list | None = None,
                 reaction_label_shift: list | None = None,
                 database: str = "BIGG",
                 remove_orphan_metabolites: bool = False,
                 
                 axis_epsilon: float = 2,):
        
        self.m_mapper = MetaNetXMapper(resource_path("metabolite_mapping.tsv"), "first")
        self.r_mapper = MetaNetXMapper(resource_path("reaction_mapping.tsv"), "first")

        self.metabolites = metabolites
        self.reactions = reactions

        self.map_main_metabolites = []
        
        self.markers_dist = markers_dist
        self.factor = scaling_factor
        self.metabolite_label_shift = (
            list(metabolite_label_shift) if metabolite_label_shift is not None else [10, 10]
        )
        self.reaction_label_shift = (
            list(reaction_label_shift) if reaction_label_shift is not None else [10, 10]
        )

        self.remove_orphan_metabolites = remove_orphan_metabolites

        # add to init
        self.h_margin = 100
        self.w_margin = 100
        self.mm_dist_part = 0.3
        self.use_const_mm_dist = False
        self.mm_dist_const = 300
        self.axis_epsilon = axis_epsilon
        self.axis_offset = 20
        self.DB = database # SEED or BIGG

        self.segments_counter = 0

        self.visible_r_indxs = []
        self.visible_m_indxs = []
        self.visible_m_ids = []
        self.m_id_to_kegg_name = {}
        self.metabolite_id2node = {}
        self.reactions_idx2kegg = {}

        self.nodes = {}
        self.text_labels = {}
        self.canvas = {
            "x": 0,
            "y": 0,
            "width": 1000,
            "height": 1000,
        }

    def build_kegg_map(self):

        escher_map = []

        description = self._generate_description("test")
        escher_map.append(description)

        model = {}

        m_desc, m2indx_dict = self._prepare_elements_descriptions(
            self.metabolites,
            self._generate_metabolite_dict,
        )

        r_desc, r2indx_dict = self._prepare_elements_descriptions(
            self.reactions,
            self._generate_reaction_dict,
        )
        r_nodes, r2node_dict = self._prepare_reactions_nodes(self.reactions)

        # prepare multimarkers between reactions and metabolites
        r_mm_nodes, r2mm_node_dict = self._prepare_reactions_multimarkers(self.reactions)

        # compose all nodes
        global_nodes_idxs = self._make_global_idxs(m2indx_dict, r2node_dict, r2mm_node_dict)
        all_nodes = self._compose_nodes(global_nodes_idxs, m_desc, r_nodes, r_mm_nodes)

        # update edges
        r_desc = self._add_edges_to_reactions_descriptions(
            self.reactions,
            r_desc,
            global_nodes_idxs,
        )
    
        model["nodes"] = {i:j for i,j in all_nodes.items() if j}
        
        model["reactions"] = {r2indx_dict[r]: r_desc[r] for r in r_desc.keys() if r2indx_dict[r]}

        model["text_labels"] = self.text_labels
        model["canvas"] = self.canvas
        
        model["nodes"], model["reactions"] = self._multiply_positions(
            model["nodes"],
            model["reactions"],
        )
        model["canvas"] = self._tune_canvas(model["nodes"], model["canvas"] )
        model["nodes"], model["reactions"] = self._align_nodes(
            model["nodes"],
            model["reactions"],
            model["canvas"],
        )

        escher_map.append(model)

        return escher_map

    def build_map(self, cobra_model):

        escher_map = []

        description = self._generate_description("test")
        escher_map.append(description)

        model = {}

        # prepare all metabolites descriptions
        m_desc, m2indx_dict = self._prepare_elements_descriptions(
            self.metabolites,
            self._generate_metabolite_dict,
        )

        # prepare all reactions descriptions with subnodes
        r_desc, r2indx_dict = self._prepare_elements_descriptions(
            self.reactions,
            self._generate_reaction_dict,
        )
        r_nodes, r2node_dict = self._prepare_reactions_nodes(self.reactions)

        # prepare multimarkers between reactions and metabolites
        r_mm_nodes, r2mm_node_dict = self._prepare_reactions_multimarkers(self.reactions)

        # compose all nodes
        global_nodes_idxs = self._make_global_idxs(m2indx_dict, r2node_dict, r2mm_node_dict)
        all_nodes = self._compose_nodes(global_nodes_idxs, m_desc, r_nodes, r_mm_nodes)

        # update edges
        r_desc = self._add_edges_to_reactions_descriptions(
            self.reactions,
            r_desc,
            global_nodes_idxs,
        )

        # extract and prepare data from model
        (
            cobra_model_metabolites,
            anti_metabolites,
            cobra_model_reactions,
            anti_reactions,
        ) = self._parse_model(cobra_model, m_desc, r_nodes)
        all_nodes, r2indx_dict = self._subtract_not_in_model_reactions(
            global_nodes_idxs,
            all_nodes,
            anti_reactions,
            r2indx_dict,
        )
        all_nodes = self._subtract_not_in_model_metabolites(
            global_nodes_idxs,
            all_nodes,
            anti_metabolites,
        )

        if self.remove_orphan_metabolites:
            all_nodes = self._remove_orphan_metabolites(all_nodes, r_desc, r2indx_dict)

        secondary_data = self._extract_secondary_metabolites(cobra_model_reactions)
        all_nodes, r_desc = self._add_secondary_metabolites(
            secondary_data,
            all_nodes,
            r_desc,
            global_nodes_idxs,
        )

        model["nodes"] = {i:j for i,j in all_nodes.items() if j}
        
        model["reactions"] = {r2indx_dict[r]: r_desc[r] for r in r_desc.keys() if r2indx_dict[r]}

        model["text_labels"] = self.text_labels
        model["canvas"] = self.canvas
        
        model["nodes"], model["reactions"] = self._multiply_positions(
            model["nodes"],
            model["reactions"],
        )
        model["canvas"] = self._tune_canvas(model["nodes"], model["canvas"] )
        model["nodes"], model["reactions"] = self._align_nodes(
            model["nodes"],
            model["reactions"],
            model["canvas"],
        )

        escher_map.append(model)

        return escher_map
    
    def _generate_description(self, name, id="default"):

        desc = {
            "map_name": name,
            "map_id": id,
            "map_description": "",
            "homepage": "https://escher.github.io",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#"
        }

        return desc
    
    def _prepare_elements_descriptions(self, elements, generation_func):

        descs = {}
        e_to_idx = {}

        for i, elem in enumerate(elements.items()):

            name, meta_data = elem

            desc = generation_func(meta_data["ids"], meta_data, name=name)
            descs[name] = desc
            e_to_idx[name] = i

        return descs, e_to_idx

    def _generate_metabolite_dict(self, ids, metabolite, name=None, primary=True):
        
        if self.DB == "BIGG":
            id = ids["BIGG"]
        elif self.DB == "SEED":
            id = ids["SEED"]
        elif self.DB == "KEGG":
            id = ids["KEGG"]

        pos = metabolite["position"]

        m_dict = {
            "node_type": "metabolite",
            "bigg_id": id,
            "name": name,
            "node_is_primary": primary,
        }

        m_dict["x"] = pos[0]
        m_dict["y"] = pos[1]
        m_dict["label_x"] = str(float(pos[0]) + self.metabolite_label_shift[0])
        m_dict["label_y"] = str(float(pos[1]) + self.metabolite_label_shift[1])

        return m_dict
    
    def _generate_reaction_dict(self, ids, reaction, name=None):

        if self.DB == "BIGG":
            id = ids["BIGG"]
        elif self.DB == "SEED":
            id = str(ids["SEED"]) + "_c0"
        elif self.DB == "KEGG":
            id = ids["KEGG"]

        reaction_dict = {
            "name": name,
            "bigg_id": id,
            "reversibility": reaction["reversibility"] in ("reversible", True),
            "gene_reaction_rule": "",
            "genes": [],
            "metabolites": [],
            "segments": {},
        }

        reaction_dict["label_x"] = str(
            float(reaction["position"][0]) + self.reaction_label_shift[0]
        )
        reaction_dict["label_y"] = str(
            float(reaction["position"][1]) + self.reaction_label_shift[1]
        )

        reaction_dict["metabolites"].extend([{"kegg_id": self.metabolites[m]["ids"]["KEGG"],
                                              "bigg_id": self.metabolites[m]["ids"]["BIGG"],
                                              "seed_id": self.metabolites[m]["ids"]["SEED"], 
                                              "coefficient": -1} 
                                             for m in reaction["substrates"].get("main", [])])
        
        reaction_dict["metabolites"].extend([{"kegg_id": self.metabolites[m]["ids"]["KEGG"],
                                              "bigg_id": self.metabolites[m]["ids"]["BIGG"],
                                              "seed_id": self.metabolites[m]["ids"]["SEED"],
                                              "coefficient": 1} 
                                             for m in reaction["products"].get("main", [])])

        return reaction_dict
    
    def _generate_node_dict(self, type, pos):

        node = {"node_type": type, "x": pos[0], "y": pos[1]}

        return node
    
    def _prepare_reactions_nodes(self, reactions):

        nodes = {}
        r_to_node = {}

        for i, rea in enumerate(reactions.items()):

            name, meta_data = rea

            node = self._generate_node_dict("midmarker", meta_data["position"])
            nodes[name] = node
            r_to_node[name] = i

        return nodes, r_to_node
    
    def _prepare_reactions_multimarkers(self, reactions):

        nodes = {}
        r_to_node = {}

        for i, rea in enumerate(reactions.items()):

            name, meta_data = rea

            in_pos, out_pos = self._calc_multimarker_positions(meta_data)

            in_node = None
            if in_pos:
                in_node = self._generate_node_dict("multimarker", in_pos)

            out_node = None
            if out_pos:
                out_node = self._generate_node_dict("multimarker", out_pos)

            nodes[name] = {"in": in_node, "out": out_node}
            r_to_node[name] = {"in": i, "out": i + len(reactions)}

        return nodes, r_to_node
    
    def _calc_mass_center(self, mets):

        mets = mets["main"]

        if mets:

            positions = [np.array(self.metabolites[m]["position"], dtype=np.float64) for m in mets]
            center = np.mean(positions, axis=0)

            return center
        
        return None
    
    # def _calc_multimarker_positions(self, reaction_data):
        
    #     reaction_pos = np.array(reaction_data["position"], dtype=np.float64)
        
    #     substrates = reaction_data.get("substrates") or []
    #     multimarker_in_pos = self._calc_multimarker_position(reaction_pos, substrates)
        
    #     products = reaction_data.get("products") or []
    #     multimarker_out_pos = self._calc_multimarker_position(reaction_pos, products)
        
    #     return multimarker_in_pos, multimarker_out_pos
    
    # def _calc_multimarker_position(self, reaction_position, metabolites):

    #     center = None
    #     if metabolites:
    #         center = self._calc_mass_center(metabolites)

    #     if center is not None:

    #         vec = center - reaction_position

    #         if self.use_const_mm_dist:
    #             norm = np.linalg.norm(vec)
    #             norm_vec = vec / norm
    #             shift = norm_vec * self.mm_dist_const

    #         else:

    #             shift = vec * self.mm_dist_part
            
    #         return (reaction_position + shift).tolist()
        
    #     return None

    def _calc_multimarker_positions(self, reaction_data):
        """Calculate both multimarker positions for a reaction."""
        
        reaction_pos = np.array(reaction_data["position"], dtype=np.float64)
    
        substrates = reaction_data.get("substrates", {})
        products = reaction_data.get("products", {})
        
        all_mets = substrates.get("main", []) + products.get("main", [])
        all_positions = [
            np.array(self.metabolites[m]["position"], dtype=np.float64)
            for m in all_mets
        ]
        
        common_axis_type, common_axis_value = self._check_metabolites_on_same_axis(all_positions)
        
        multimarker_in_pos = self._calc_multimarker_position(
            reaction_pos, substrates, products, common_axis_type, common_axis_value
        )
        
        multimarker_out_pos = self._calc_multimarker_position(
            reaction_pos, products, substrates, common_axis_type, common_axis_value
        )
        
        return multimarker_in_pos, multimarker_out_pos
    
    def _calc_multimarker_position(self, reaction_pos, metabolites, opposite_metabolites, 
                                common_axis_type=None, common_axis_value=None):
        """
        Calculate one multimarker position.
        """
        
        mets = metabolites.get("main", None)
        if not mets:
            return None
        
        positions = [np.array(self.metabolites[m]["position"], dtype=np.float64) for m in mets]
        
        if common_axis_type == "vertical":
            mm_x = reaction_pos[0]
            my_center_y = np.mean([p[1] for p in positions])
            offset_dir = np.sign(my_center_y - reaction_pos[1])
            offset_dir = offset_dir if offset_dir != 0 else 1
            mm_y = reaction_pos[1] + self.axis_offset * offset_dir
            return [mm_x, mm_y]
        
        if common_axis_type == "horizontal":
            mm_y = reaction_pos[1]
            my_center_x = np.mean([p[0] for p in positions])
            offset_dir = np.sign(my_center_x - reaction_pos[0])
            offset_dir = offset_dir if offset_dir != 0 else 1
            mm_x = reaction_pos[0] + self.axis_offset * offset_dir
            return [mm_x, mm_y]
        
        aligned_type, aligned_pos = self._find_aligned_metabolite(positions, reaction_pos)
        
        if aligned_type is not None:
            # TODO: Revisit opposite-side context during the layout refactor.
            opposite_mets = opposite_metabolites.get("main", [])  # noqa: F841
            
            if aligned_type == "horizontal":
                mm_y = reaction_pos[1]
                
                offset_dir = np.sign(aligned_pos[0] - reaction_pos[0])
                offset_dir = offset_dir if offset_dir != 0 else 1
                mm_x = reaction_pos[0] + self.axis_offset * offset_dir
                
                return [mm_x, mm_y]

            if aligned_type == "vertical":
                mm_x = reaction_pos[0]
                
                offset_dir = np.sign(aligned_pos[1] - reaction_pos[1])
                offset_dir = offset_dir if offset_dir != 0 else 1
                mm_y = reaction_pos[1] + self.axis_offset * offset_dir
                
                return [mm_x, mm_y]
        
        return self._calc_multimarker_position_by_mass_center(reaction_pos, positions)
    
    def _check_metabolites_on_same_axis(self, positions):
        """
        Return the shared axis for aligned metabolite positions, if any.
        """
        if len(positions) < 2:
            return None, None
        
        ys = [p[1] for p in positions]
        xs = [p[0] for p in positions]
        
        if max(ys) - min(ys) < self.axis_epsilon:
            return "horizontal", np.mean(ys)
        
        if max(xs) - min(xs) < self.axis_epsilon:
            return "vertical", np.mean(xs)
        
        return None, None
    
    def _find_aligned_metabolite(self, positions, reaction_pos):
        best_dist = float('inf')
        best_type = None
        best_pos = None
        
        for pos in positions:
            dx = abs(pos[0] - reaction_pos[0])
            dy = abs(pos[1] - reaction_pos[1])
            dist = np.linalg.norm(pos - reaction_pos)
            
            on_horizontal = dy < self.axis_epsilon
            on_vertical = dx < self.axis_epsilon
            
            if on_horizontal and on_vertical:
                if dy < dx:
                    align_type = "horizontal"
                else:
                    align_type = "vertical"
            elif on_horizontal:
                align_type = "horizontal"
            elif on_vertical:
                align_type = "vertical"
            else:
                continue
            
            if dist < best_dist:
                best_dist = dist
                best_type = align_type
                best_pos = pos
        
        return best_type, best_pos
    
    def _find_nearest_opposite_coord(self, reaction_pos, opposite_mets, coord_idx):
        """
        Return the selected coordinate from the nearest opposite metabolite.
        """
        best_dist = float('inf')
        best_coord = reaction_pos[coord_idx]
        
        for m in opposite_mets:
            pos = np.array(self.metabolites[m]["position"], dtype=np.float64)
            dist = np.linalg.norm(pos - reaction_pos)
            
            if dist < best_dist:
                best_dist = dist
                best_coord = pos[coord_idx]
        
        return best_coord
    
    def _calc_multimarker_position_by_mass_center(self, reaction_pos, positions):
        """Calculate fallback position from the mass center."""
        
        center = np.mean(positions, axis=0)
        vec = center - reaction_pos
        
        if self.use_const_mm_dist:
            norm = np.linalg.norm(vec)
            if norm > 0:
                norm_vec = vec / norm
                shift = norm_vec * self.mm_dist_const
            else:
                shift = np.array([self.mm_dist_const, 0])
        else:
            shift = vec * self.mm_dist_part
        
        return (reaction_pos + shift).tolist()

    def _make_global_idxs(self, m2i, r2i, mm2i):

        global_idxs = {}

        global_idxs["metabolites"] = {m: m2i[m] for m in m2i.keys()}
        global_idxs["reactions"] = {r: r2i[r] + len(m2i) for r in r2i.keys()}

        mr_len = len(m2i) + len(r2i)

        global_idxs["multimarkers"] = {r: {"in": mr_len + mm2i[r]["in"],
                                            "out": mr_len + mm2i[r]["out"]} 
                                       for r in mm2i.keys()}
    
        return global_idxs
    
    def _compose_nodes(self, global_idxs, m_nodes, r_nodes, mm_nodes):

        nodes = {}

        for m, node in m_nodes.items():

            nodes[global_idxs["metabolites"][m]] = node

        for r, node in r_nodes.items():

            nodes[global_idxs["reactions"][r]] = node

        for r, in_out_nodes in mm_nodes.items():

            if in_out_nodes["in"]:
                nodes[global_idxs["multimarkers"][r]["in"]] = in_out_nodes["in"]
            if in_out_nodes["out"]:
                nodes[global_idxs["multimarkers"][r]["out"]] = in_out_nodes["out"]

        return nodes
    
    def _add_edges_to_reactions_descriptions(self, reactions, reactions_descs, global_idxs):

        edges_num = 0

        for r, meta_data in reactions.items():

            r_edges = self._prepare_reaction_edges(r, meta_data, global_idxs)
            reactions_descs[r]["segments"] = {edges_num + i: edge for i, edge in enumerate(r_edges)}

            edges_num += len(r_edges)

        return reactions_descs

    def _prepare_reaction_edges(self, r_name, reaction, global_idxs):

        edges = []

        in_mm_node_idx = global_idxs["multimarkers"][r_name]["in"]
        out_mm_node_idx = global_idxs["multimarkers"][r_name]["out"]
        reaction_node_idx = global_idxs["reactions"][r_name]

        for m in reaction["substrates"].get("main", []):
            edges.append(self._prepare_edge_dict(global_idxs["metabolites"][m], in_mm_node_idx))

        edges.append(self._prepare_edge_dict(in_mm_node_idx, reaction_node_idx))
        edges.append(self._prepare_edge_dict(reaction_node_idx, out_mm_node_idx))

        for m in reaction["products"].get("main", []):
            edges.append(self._prepare_edge_dict(out_mm_node_idx, global_idxs["metabolites"][m]))

        return edges

    def _prepare_edge_dict(self, from_node, to_node):

        edge_dict = {
            "from_node_id": from_node,
            "to_node_id": to_node,
            "b1": None,
            "b2": None,
        }

        return edge_dict
    
    # model integration

    def _parse_model(self, model, m_nodes, r_nodes):

        matched_rs, anti_rs = self._extract_model_reactions(model, r_nodes)
        ms, anti_ms = self._extract_model_metabolites(model, m_nodes)

        return ms, anti_ms, matched_rs, anti_rs
    
    def _extract_model_reactions(self, model, r_nodes):

        matched = {}
        anti_reactions = []

        for r_name in r_nodes.keys():

            keggs = set([r_name])
            biggs = set(self.r_mapper[r_name].bigg_all) if self.r_mapper.get(r_name) else set()
            seeds = set(self.r_mapper[r_name].seed_all) if self.r_mapper.get(r_name) else set()

            found = None
            for rxn in model.reactions:

                rxn_kegg = rxn.annotation.get("kegg.reaction", [])
                rxn_kegg = rxn_kegg if isinstance(rxn_kegg, list) else [rxn_kegg]

                rxn_bigg = rxn.annotation.get("bigg.reaction", [])
                rxn_bigg = rxn_bigg if isinstance(rxn_bigg, list) else [rxn_bigg]

                rxn_seed = rxn.annotation.get("seed.reaction", [])
                rxn_seed = rxn_seed if isinstance(rxn_seed, list) else [rxn_seed]

                if (keggs & set(rxn_kegg) or
                    biggs & set(rxn_bigg) or
                    seeds & set(rxn_seed)):
                    found = rxn
                    break

            if found:
                matched[r_name] = found
            else:
                anti_reactions.append(r_name)

        return matched, anti_reactions

    def _extract_model_metabolites(self, model, m_nodes):

        model_ids = {"KEGG": set(), "BIGG": set(), "SEED": set()}

        for met in model.metabolites:
            for k in model_ids.keys():
                key = f'{k.lower()}.{"metabolite" if k == "BIGG" else "compound"}'
                names = met.annotation.get(key, [])
                names = names if isinstance(names, list) else [names]
                model_ids[k] |= set(names)

        ms = {}
        anti_ms = []

        for m_name in m_nodes.keys():

            keggs = {m_name}
            biggs = set(self.m_mapper[m_name].bigg_all) if self.m_mapper.get(m_name) else set()
            seeds = set(self.m_mapper[m_name].seed_all) if self.m_mapper.get(m_name) else set()

            if keggs & model_ids["KEGG"] or biggs & model_ids["BIGG"] or seeds & model_ids["SEED"]:
                ms[m_name] = True
            else:
                anti_ms.append(m_name)

        return ms, anti_ms
    
    def _subtract_not_in_model_reactions(self, global_idxs, all_nodes, anti_rs, r2indx_dict):

        for r_name in anti_rs:

            all_nodes[global_idxs["reactions"][r_name]] = None
            all_nodes[global_idxs["multimarkers"][r_name]["in"]] = None
            all_nodes[global_idxs["multimarkers"][r_name]["out"]] = None

            r2indx_dict[r_name] = None

        return all_nodes, r2indx_dict
    
    def _subtract_not_in_model_metabolites(self, global_idxs, nodes, anti_ms):

        for m in anti_ms:

            nodes[global_idxs["metabolites"][m]] = None

        return nodes
    
    def _extract_secondary_metabolites(self, matched_rs):
        """
        Extract model metabolites that are not primary map metabolites.
        
        Returns {r_name: {"substrates": [...], "products": [...]}}.
        """

        main_met_ids = set()
        for _m_name, m_data in self.metabolites.items():
            ids = m_data["ids"]
            for v in ids.values():
                if v:
                    main_met_ids.add(v)

        secondary = {}

        for r_name, cobra_rxn in matched_rs.items():

            sec_subs = []
            sec_prods = []

            for met, coef in cobra_rxn.metabolites.items():

                met_ids = set()
                met_ids.add(met.id)
                met_ids.add(met.id[:-2] if len(met.id) > 2 else met.id)

                for key in ["kegg.compound", "bigg.metabolite", "seed.compound"]:
                    ann = met.annotation.get(key, [])
                    if isinstance(ann, str):
                        ann = [ann]
                    met_ids.update(ann)

                if met_ids & main_met_ids:
                    continue

                entry = {"bigg_id": met.id, "name": met.name, "coefficient": coef}

                if coef < 0:
                    sec_subs.append(entry)
                else:
                    sec_prods.append(entry)

            if sec_subs or sec_prods:
                secondary[r_name] = {
                    "substrates": sec_subs,
                    "products": sec_prods,
                }

        return secondary
    
    def _add_secondary_metabolites(self, secondary_data, all_nodes, r_desc, global_idxs):

        max_node_idx = max(int(k) for k in all_nodes.keys()) + 1
        
        seg_counter = 0
        for _r_name, r_data in r_desc.items():
            if r_data["segments"]:
                seg_counter = max(seg_counter, max(int(k) for k in r_data["segments"].keys()) + 1)

        for r_name, sec in secondary_data.items():

            r_node_idx = global_idxs["reactions"][r_name]
            in_mm_idx = global_idxs["multimarkers"][r_name]["in"]
            out_mm_idx = global_idxs["multimarkers"][r_name]["out"]

            if (all_nodes.get(r_node_idx) is None or 
                all_nodes.get(in_mm_idx) is None or 
                all_nodes.get(out_mm_idx) is None):
                continue

            reaction_pos = np.array(
                [all_nodes[r_node_idx]["x"], all_nodes[r_node_idx]["y"]],
                dtype=np.float64,
            )
            in_mm_pos = np.array(
                [all_nodes[in_mm_idx]["x"], all_nodes[in_mm_idx]["y"]],
                dtype=np.float64,
            )
            out_mm_pos = np.array(
                [all_nodes[out_mm_idx]["x"], all_nodes[out_mm_idx]["y"]],
                dtype=np.float64,
            )

            subs_center = self._calc_main_metabolites_center(
                self.reactions[r_name]["substrates"].get("main", [])
            )
            prods_center = self._calc_main_metabolites_center(
                self.reactions[r_name]["products"].get("main", [])
            )

            if sec["substrates"] and subs_center is not None:
                direction, perp = self._calc_secondary_directions_from_center(
                    reaction_pos,
                    subs_center,
                )
                for j, entry in enumerate(sec["substrates"]):
                    pos = self._calc_secondary_position(
                        in_mm_pos,
                        direction,
                        perp,
                        j,
                        len(sec["substrates"]),
                        side=1,
                    )
                    node = self._generate_secondary_metabolite_dict(
                        entry["bigg_id"],
                        entry["name"],
                        pos,
                    )
                    all_nodes[max_node_idx] = node
                    r_desc[r_name]["segments"][seg_counter] = self._prepare_edge_dict(
                        max_node_idx,
                        in_mm_idx,
                    )
                    r_desc[r_name]["metabolites"].append(
                        {
                            "bigg_id": entry["bigg_id"],
                            "coefficient": entry["coefficient"],
                        }
                    )
                    seg_counter += 1
                    max_node_idx += 1

            if sec["products"] and prods_center is not None:
                direction, perp = self._calc_secondary_directions_from_center(
                    reaction_pos,
                    prods_center,
                )
                for j, entry in enumerate(sec["products"]):
                    pos = self._calc_secondary_position(
                        out_mm_pos,
                        direction,
                        perp,
                        j,
                        len(sec["products"]),
                        side=1,
                    )
                    node = self._generate_secondary_metabolite_dict(
                        entry["bigg_id"],
                        entry["name"],
                        pos,
                    )
                    all_nodes[max_node_idx] = node
                    r_desc[r_name]["segments"][seg_counter] = self._prepare_edge_dict(
                        out_mm_idx,
                        max_node_idx,
                    )
                    r_desc[r_name]["metabolites"].append(
                        {
                            "bigg_id": entry["bigg_id"],
                            "coefficient": entry["coefficient"],
                        }
                    )
                    seg_counter += 1
                    max_node_idx += 1

        return all_nodes, r_desc


    def _calc_main_metabolites_center(self, met_names):
        
        if not met_names:
            return None
        
        positions = [np.array(self.metabolites[m]["position"], dtype=np.float64) for m in met_names]
        return np.mean(positions, axis=0)


    def _calc_secondary_directions_from_center(self, reaction_pos, mets_center):
        
        vec = mets_center - reaction_pos
        norm = np.linalg.norm(vec)
        
        if norm > 0:
            direction = vec / norm
        else:
            direction = np.array([1.0, 0.0])
        
        perp = np.array([-direction[1], direction[0]])
        
        return direction, perp
        
    def _calc_secondary_directions(self, in_mm_pos, out_mm_pos):

        reaction_vec = out_mm_pos - in_mm_pos
        norm = np.linalg.norm(reaction_vec)
        
        if norm > 0:
            reaction_dir = reaction_vec / norm
        else:
            reaction_dir = np.array([1.0, 0.0])

        perp = np.array([-reaction_dir[1], reaction_dir[0]])

        return reaction_dir, perp


    def _calc_secondary_position(self, anchor_pos, reaction_dir, perp, index, total, side=-1):
        """
        side: -1 for substrates, +1 for products.
        """

        lateral_offset = (-(total - 1) / 2.0 + index) * self.markers_dist * 3
        pos = anchor_pos + side * reaction_dir * self.markers_dist * 2 + perp * lateral_offset

        return pos.tolist()


    def _generate_secondary_metabolite_dict(self, bigg_id, name, pos):

        return {
            "node_type": "metabolite",
            "bigg_id": bigg_id,
            "name": name,
            "node_is_primary": False,
            "x": pos[0],
            "y": pos[1],
            "label_x": pos[0] + self.metabolite_label_shift[0],
            "label_y": pos[1] + self.metabolite_label_shift[1],
        }
    
    def _remove_orphan_metabolites(self, all_nodes, r_desc, r2indx_dict):

        referenced_nodes = set()
        for r_name, r_data in r_desc.items():
            if r2indx_dict.get(r_name) is None:
                continue
            for seg in r_data["segments"].values():
                referenced_nodes.add(seg["from_node_id"])
                referenced_nodes.add(seg["to_node_id"])

        for nid, node in all_nodes.items():
            if (
                node
                and node.get("node_type") == "metabolite"
                and node.get("node_is_primary")
                and nid not in referenced_nodes
            ):
                all_nodes[nid] = None

        return all_nodes

    # scaling and canvas 

    def _tune_canvas(self, nodes, canvas):

        x, y = 0, 0
        min_x, min_y = 0, 0

        for _i, compound_dict in nodes.items():

            if not compound_dict:
                continue

            if compound_dict["x"]:

                if x < float(compound_dict["x"]):
                    x = float(compound_dict["x"])
                if y < float(compound_dict["y"]):
                    y = float(compound_dict["y"])

                if min_x > float(compound_dict["x"]):
                    min_x = float(compound_dict["x"])
                if min_y > float(compound_dict["y"]):
                    min_y = float(compound_dict["y"])

        canvas["width"] = x - min_x + self.h_margin
        canvas["height"] = y - min_y + self.w_margin

        return canvas

    def _multiply_positions(self, nodes, reactions):

        for node in nodes.values():

            if not node:
                continue

            node["x"] = float(node["x"]) * self.factor
            node["y"] = float(node["y"]) * self.factor

            if "label_x" in node.keys():

                node["label_x"] = float(node["label_x"]) * self.factor
                node["label_y"] = float(node["label_y"]) * self.factor

        for reaction in reactions.values():

            if not reaction:
                continue

            reaction["label_x"] = float(reaction["label_x"]) * self.factor
            reaction["label_y"] = float(reaction["label_y"]) * self.factor

        return nodes, reactions
    
    def _align_nodes(self, nodes, reactions, canvas):

        canvas_x, canvas_y = canvas["width"] / 2, canvas["height"] / 2
        current_x, current_y = self._current_center(nodes, reactions)

        shift_x, shift_y = canvas_x - current_x, canvas_y - current_y

        for node in nodes.values():

            if not node:
                continue

            node["x"] = float(node["x"]) + shift_x
            node["y"] = float(node["y"]) + shift_y

            if "label_x" in node.keys():

                node["label_x"] = float(node["label_x"]) + shift_x
                node["label_y"] = float(node["label_y"]) + shift_y

        for reaction in reactions.values():

            if not reaction:
                continue

            reaction["label_x"] = float(reaction["label_x"]) + shift_x
            reaction["label_y"] = float(reaction["label_y"]) + shift_y

        return nodes, reactions

    def _current_center(self, nodes, reactions):

        min_x, min_y = float('inf'), float('inf')
        max_x, max_y = float('-inf'), float('-inf')

        for node in nodes.values():
            if not node:
                continue
            min_x = min(min_x, float(node["x"]))
            max_x = max(max_x, float(node["x"]))
            min_y = min(min_y, float(node["y"]))
            max_y = max(max_y, float(node["y"]))

        return (min_x + max_x) / 2, (min_y + max_y) / 2
    
