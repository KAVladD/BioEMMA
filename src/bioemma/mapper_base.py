import numpy as np

from bioemma.metanetx_mapper import MetaNetXMapper
from bioemma._resources import resource_path

class EscherMapper:

    def __init__(self, 
                 metabolites: dict, 
                 reactions: dict,
                 markers_dist: float = 10,
                 scaling_factor: float = 4,
                 metabolite_label_shift: list | None = None,
                 reaction_label_shift: list | None = None,
                 database: str = "BIGG",
                 remove_orphan_metabolites: bool = False,
                 include_kegg_only: bool = False,
                 canvas_margin_x: float = 160,
                 canvas_margin_y: float = 160,
                 multimarker_distance_fraction: float = 0.3,
                 use_constant_multimarker_distance: bool = False,
                 constant_multimarker_distance: float = 300,
                 axis_offset: float = 20,
                 secondary_metabolite_distance: float | None = None,
                 secondary_metabolite_spacing: float | None = None,
                 use_model_metabolite_ids: bool = False,
                 use_database_secondary_metabolite_ids: bool = False,
                 metabolite_id_compartments: bool | None = None,
                 compartment_filter: str | None = None,
                 canvas_width: float = 1000,
                 canvas_height: float = 1000,
                 axis_epsilon: float = 2,):
        
        self.m_mapper = MetaNetXMapper(resource_path("metabolite_mapping.tsv"), "first")
        self.r_mapper = MetaNetXMapper(resource_path("reaction_mapping.tsv"), "first")

        self.metabolites = metabolites
        self.reactions = reactions

        self.map_main_metabolites = []
        
        self.markers_dist = float(markers_dist)
        self.factor = scaling_factor
        self.metabolite_label_shift = (
            list(metabolite_label_shift) if metabolite_label_shift is not None else [10, 10]
        )
        self.reaction_label_shift = (
            list(reaction_label_shift) if reaction_label_shift is not None else [10, 10]
        )

        self.remove_orphan_metabolites = remove_orphan_metabolites
        self.include_kegg_only = include_kegg_only

        self.canvas_margin_x = float(canvas_margin_x)
        self.canvas_margin_y = float(canvas_margin_y)
        self.h_margin = self.canvas_margin_x
        self.w_margin = self.canvas_margin_y
        self.mm_dist_part = float(multimarker_distance_fraction)
        self.use_const_mm_dist = use_constant_multimarker_distance
        self.mm_dist_const = float(constant_multimarker_distance)
        self.axis_epsilon = axis_epsilon
        self.axis_offset = float(axis_offset)
        self.secondary_metabolite_distance = (
            float(secondary_metabolite_distance)
            if secondary_metabolite_distance is not None
            else self.markers_dist * 2
        )
        self.secondary_metabolite_spacing = (
            float(secondary_metabolite_spacing)
            if secondary_metabolite_spacing is not None
            else self.markers_dist * 3
        )
        self.DB = database # SEED or BIGG
        self.use_model_metabolite_ids = use_model_metabolite_ids
        self.use_database_secondary_metabolite_ids = (
            use_database_secondary_metabolite_ids
        )
        self.metabolite_id_compartments = (
            self.DB == "BIGG"
            if metabolite_id_compartments is None
            else bool(metabolite_id_compartments)
        )
        self.compartment_filter = self._normalize_compartment_filter(
            compartment_filter
        )

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
            "width": canvas_width,
            "height": canvas_height,
        }
        self.map_stats = {}

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
        
        model["reactions"] = {
            r2indx_dict[r]: r_desc[r]
            for r in r_desc.keys()
            if r2indx_dict[r] is not None
        }

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

        self.map_stats = {"stages": []}
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
        self._record_map_stage(
            "kegg_layout",
            all_nodes,
            r_desc,
            r2indx_dict,
            description="Initial Escher layout reconstructed from KEGG.",
        )

        # extract and prepare data from model
        (
            cobra_model_metabolites,
            anti_metabolites,
            cobra_model_reactions,
            anti_reactions,
            model_reaction_bigg_ids,
        ) = self._parse_model(cobra_model, m_desc, r_nodes)
        self._apply_model_reaction_bigg_ids(r_desc, model_reaction_bigg_ids)
        if self.use_model_metabolite_ids:
            self._apply_model_metabolite_ids(
                all_nodes,
                r_desc,
                global_nodes_idxs,
                cobra_model_metabolites,
            )
            self._record_map_stage(
                "model_metabolite_id_application",
                all_nodes,
                r_desc,
                r2indx_dict,
                description=(
                    "Applied COBRA model metabolite identifiers to matched "
                    "KEGG metabolites."
                ),
            )
        self.map_stats["model_matching"] = {
            "matched_metabolites": len(cobra_model_metabolites),
            "unmatched_metabolites": len(anti_metabolites),
            "matched_reactions": len(cobra_model_reactions),
            "unmatched_reactions": len(anti_reactions),
            "compartment_filter": self.compartment_filter,
        }
        if not self.include_kegg_only:
            all_nodes, r2indx_dict = self._subtract_not_in_model_reactions(
                global_nodes_idxs,
                all_nodes,
                anti_reactions,
                r2indx_dict,
            )
            self._record_map_stage(
                "model_reaction_filter",
                all_nodes,
                r_desc,
                r2indx_dict,
                description="Removed KEGG reactions that were not matched to the COBRA model.",
            )
            all_nodes = self._subtract_not_in_model_metabolites(
                global_nodes_idxs,
                all_nodes,
                anti_metabolites,
                r_desc,
                r2indx_dict,
            )
            self._record_map_stage(
                "model_metabolite_filter",
                all_nodes,
                r_desc,
                r2indx_dict,
                description="Removed KEGG metabolites that were not matched to the COBRA model.",
            )
        else:
            self._record_map_stage(
                "model_filter_skipped",
                all_nodes,
                r_desc,
                r2indx_dict,
                description="Kept KEGG-only reactions and metabolites.",
            )

        secondary_data = self._extract_secondary_metabolites(cobra_model_reactions)
        all_nodes, r_desc = self._add_secondary_metabolites(
            secondary_data,
            all_nodes,
            r_desc,
            global_nodes_idxs,
        )
        self._record_map_stage(
            "secondary_metabolite_addition",
            all_nodes,
            r_desc,
            r2indx_dict,
            description="Added non-primary COBRA metabolites attached to matched reactions.",
        )

        if self.remove_orphan_metabolites:
            all_nodes = self._remove_orphan_metabolites(all_nodes, r_desc, r2indx_dict)
            self._record_map_stage(
                "orphan_metabolite_filter",
                all_nodes,
                r_desc,
                r2indx_dict,
                description=(
                    "Removed metabolite nodes that are not referenced "
                    "by any visible reaction."
                ),
            )

        model["nodes"] = {i:j for i,j in all_nodes.items() if j}
        
        model["reactions"] = {
            r2indx_dict[r]: r_desc[r]
            for r in r_desc.keys()
            if r2indx_dict[r] is not None
        }

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
        self._record_map_stage(
            "final_layout",
            model["nodes"],
            model["reactions"],
            description="Final Escher map after layout scaling and alignment.",
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

    def _record_map_stage(
        self,
        name,
        nodes,
        reactions,
        reaction_index=None,
        description=None,
    ):
        counts = self._count_map_elements(nodes, reactions, reaction_index)
        previous_counts = (
            self.map_stats["stages"][-1]["counts"]
            if self.map_stats.get("stages")
            else None
        )
        entry = {
            "name": name,
            "counts": counts,
            "change": self._count_delta(previous_counts, counts),
        }
        if description:
            entry["description"] = description
        self.map_stats.setdefault("stages", []).append(entry)

    def _count_map_elements(self, nodes, reactions, reaction_index=None):
        visible_nodes = [node for node in nodes.values() if node]
        if reaction_index is None:
            visible_reactions = [reaction for reaction in reactions.values() if reaction]
        else:
            visible_reactions = [
                reactions[name]
                for name in reactions
                if reaction_index.get(name) is not None and reactions[name]
            ]

        segments = sum(len(reaction.get("segments", {})) for reaction in visible_reactions)
        primary_metabolites = sum(
            1
            for node in visible_nodes
            if node.get("node_type") == "metabolite" and node.get("node_is_primary")
        )
        secondary_metabolites = sum(
            1
            for node in visible_nodes
            if node.get("node_type") == "metabolite"
            and node.get("node_is_primary") is False
        )
        reaction_nodes = sum(1 for node in visible_nodes if node.get("node_type") == "midmarker")
        multimarkers = sum(1 for node in visible_nodes if node.get("node_type") == "multimarker")

        return {
            "total_elements": len(visible_nodes) + len(visible_reactions) + segments,
            "nodes": len(visible_nodes),
            "primary_metabolites": primary_metabolites,
            "secondary_metabolites": secondary_metabolites,
            "reaction_nodes": reaction_nodes,
            "multimarkers": multimarkers,
            "reactions": len(visible_reactions),
            "segments": segments,
        }

    def _count_delta(self, before, after):
        if before is None:
            return {
                key: {"added": value, "removed": 0, "delta": value}
                for key, value in after.items()
            }

        delta = {}
        for key, value in after.items():
            diff = value - before.get(key, 0)
            delta[key] = {
                "added": max(diff, 0),
                "removed": max(-diff, 0),
                "delta": diff,
            }
        return delta
    
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
            offset = self._calc_multimarker_distance(
                abs(my_center_y - reaction_pos[1])
            )
            mm_y = reaction_pos[1] + offset * offset_dir
            return [mm_x, mm_y]
        
        if common_axis_type == "horizontal":
            mm_y = reaction_pos[1]
            my_center_x = np.mean([p[0] for p in positions])
            offset_dir = np.sign(my_center_x - reaction_pos[0])
            offset_dir = offset_dir if offset_dir != 0 else 1
            offset = self._calc_multimarker_distance(
                abs(my_center_x - reaction_pos[0])
            )
            mm_x = reaction_pos[0] + offset * offset_dir
            return [mm_x, mm_y]
        
        aligned_type, aligned_pos = self._find_aligned_metabolite(positions, reaction_pos)
        
        if aligned_type is not None:
            # TODO: Revisit opposite-side context during the layout refactor.
            opposite_mets = opposite_metabolites.get("main", [])  # noqa: F841
            
            if aligned_type == "horizontal":
                mm_y = reaction_pos[1]
                
                offset_dir = np.sign(aligned_pos[0] - reaction_pos[0])
                offset_dir = offset_dir if offset_dir != 0 else 1
                offset = self._calc_multimarker_distance(
                    abs(aligned_pos[0] - reaction_pos[0])
                )
                mm_x = reaction_pos[0] + offset * offset_dir
                
                return [mm_x, mm_y]

            if aligned_type == "vertical":
                mm_x = reaction_pos[0]
                
                offset_dir = np.sign(aligned_pos[1] - reaction_pos[1])
                offset_dir = offset_dir if offset_dir != 0 else 1
                offset = self._calc_multimarker_distance(
                    abs(aligned_pos[1] - reaction_pos[1])
                )
                mm_y = reaction_pos[1] + offset * offset_dir
                
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

    def _calc_multimarker_distance(self, base_distance):
        if self.use_const_mm_dist:
            return self.mm_dist_const

        offset = base_distance * self.mm_dist_part
        if offset > 0:
            return offset

        return self.axis_offset
    
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
        
        norm = np.linalg.norm(vec)
        if norm > 0:
            norm_vec = vec / norm
            shift = norm_vec * self._calc_multimarker_distance(norm)
        else:
            shift = np.array([self._calc_multimarker_distance(0), 0])
        
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

        matched_rs, anti_rs, model_reaction_bigg_ids = self._extract_model_reactions(
            model,
            r_nodes,
        )
        ms, anti_ms = self._extract_model_metabolites(model, m_nodes, matched_rs)

        return ms, anti_ms, matched_rs, anti_rs, model_reaction_bigg_ids
    
    def _extract_model_reactions(self, model, r_nodes):

        matched = {}
        model_reaction_bigg_ids = {}
        anti_reactions = []

        for r_name in r_nodes.keys():

            keggs = set([r_name])
            mapping_entry = self.r_mapper.get(r_name)
            biggs = set(mapping_entry.bigg_all) if mapping_entry else set()
            seeds = set(mapping_entry.seed_all) if mapping_entry else set()
            metacycs = set(mapping_entry.metacyc_all) if mapping_entry else set()
            rheas = set(mapping_entry.rhea_all) if mapping_entry else set()
            ecs = set(mapping_entry.ec_all) if mapping_entry else set()

            best_match = None
            best_score = -1
            best_bigg_id = None
            for rxn in model.reactions:
                if not self._reaction_matches_compartment_filter(rxn):
                    continue

                rxn_kegg = self._annotation_values(rxn.annotation, "kegg.reaction")
                rxn_bigg = self._annotation_values(rxn.annotation, "bigg.reaction")
                rxn_seed = self._annotation_values(rxn.annotation, "seed.reaction")
                rxn_metacyc = self._annotation_values(
                    rxn.annotation,
                    "metacyc.reaction",
                )
                rxn_rhea = self._annotation_values(rxn.annotation, "rhea")
                rxn_ec = self._annotation_values_any(
                    rxn.annotation,
                    ("ec-code", "ec.number", "ec_number", "ec"),
                )

                bigg_matches = biggs & set(rxn_bigg)
                score = 0
                candidate_bigg_id = None

                if bigg_matches:
                    score = 100
                    candidate_bigg_id = self._select_model_reaction_bigg_id(
                        rxn_bigg,
                        bigg_matches,
                    )
                elif keggs & set(rxn_kegg):
                    score = 90
                elif seeds & set(rxn_seed):
                    score = 80
                elif metacycs & set(rxn_metacyc):
                    score = 70
                elif rheas & set(rxn_rhea):
                    score = 70
                elif ecs & set(rxn_ec):
                    score = 10

                if score and score > best_score:
                    best_match = rxn
                    best_score = score
                    best_bigg_id = candidate_bigg_id

            if best_match:
                matched[r_name] = best_match
                if best_bigg_id:
                    model_reaction_bigg_ids[r_name] = best_bigg_id
            else:
                anti_reactions.append(r_name)

        return matched, anti_reactions, model_reaction_bigg_ids

    def _annotation_values(self, annotation, key):
        values = annotation.get(key, [])
        if values is None:
            return []
        if isinstance(values, str):
            values = [values]
        elif not isinstance(values, (list, tuple, set)):
            values = [values]
        return [str(value) for value in values if value]

    def _annotation_values_any(self, annotation, keys):
        values = []
        for key in keys:
            values.extend(self._annotation_values(annotation, key))
        return values

    def _select_model_reaction_bigg_id(self, rxn_bigg, bigg_matches):
        for value in rxn_bigg:
            if value in bigg_matches:
                return value
        return sorted(bigg_matches)[0] if bigg_matches else None

    def _apply_model_reaction_bigg_ids(self, r_desc, model_reaction_bigg_ids):
        for r_name, bigg_id in model_reaction_bigg_ids.items():
            if r_name in r_desc:
                r_desc[r_name]["bigg_id"] = bigg_id

    def _apply_model_metabolite_ids(
        self,
        all_nodes,
        r_desc,
        global_idxs,
        model_metabolites,
    ):
        for m_name, model_metabolite in model_metabolites.items():
            model_id = self._model_metabolite_output_id(model_metabolite, m_name)
            compartment = self._model_metabolite_compartment(model_metabolite)

            node_idx = global_idxs["metabolites"].get(m_name)
            if node_idx is not None and all_nodes.get(node_idx) is not None:
                all_nodes[node_idx]["bigg_id"] = model_id
                if self.metabolite_id_compartments and compartment:
                    all_nodes[node_idx]["compartment"] = compartment

            map_ids = self._map_metabolite_ids(m_name)
            for reaction in r_desc.values():
                for reaction_metabolite in reaction.get("metabolites", []):
                    values = {
                        reaction_metabolite.get("kegg_id"),
                        reaction_metabolite.get("bigg_id"),
                        reaction_metabolite.get("seed_id"),
                    }
                    if values & map_ids:
                        reaction_metabolite["bigg_id"] = model_id
                        if self.metabolite_id_compartments and compartment:
                            reaction_metabolite["compartment"] = compartment

    def _extract_model_metabolites(self, model, m_nodes, matched_rs=None):

        ms = {}
        anti_ms = []

        for m_name in m_nodes.keys():
            candidates = self._find_connected_model_metabolite_candidates(
                m_name,
                matched_rs or {},
            )
            if not candidates:
                candidates = [
                    (met, 1)
                    for met in model.metabolites
                    if self._model_metabolite_matches_map_metabolite(met, m_name)
                    and self._model_metabolite_matches_compartment_filter(met)
                ]

            if candidates:
                ms[m_name] = self._select_model_metabolite_candidate(candidates)
            else:
                anti_ms.append(m_name)

        return ms, anti_ms

    def _find_connected_model_metabolite_candidates(self, m_name, matched_rs):
        candidate_counts = {}

        for r_name, cobra_rxn in matched_rs.items():
            reaction = self.reactions.get(r_name, {})
            primary_metabolites = set(reaction.get("substrates", {}).get("main", []))
            primary_metabolites.update(reaction.get("products", {}).get("main", []))
            if m_name not in primary_metabolites:
                continue

            for met in cobra_rxn.metabolites:
                if self._model_metabolite_matches_map_metabolite(met, m_name):
                    candidate_counts[met] = candidate_counts.get(met, 0) + 1

        return list(candidate_counts.items())

    def _select_model_metabolite_candidate(self, candidates):
        return sorted(
            candidates,
            key=lambda item: (
                -item[1],
                not self._model_metabolite_matches_compartment_filter(item[0]),
                item[0].id,
            ),
        )[0][0]

    def _model_metabolite_matches_map_metabolite(self, met, m_name):
        return bool(self._model_metabolite_ids(met) & self._map_metabolite_ids(m_name))

    def _map_metabolite_ids(self, m_name):
        ids = {m_name}
        m_data = self.metabolites.get(m_name, {})
        for value in m_data.get("ids", {}).values():
            if value:
                ids.add(str(value))

        mapped = self.m_mapper.get(m_name)
        if mapped:
            ids.update(str(value) for value in mapped.bigg_all if value)
            ids.update(str(value) for value in mapped.seed_all if value)

        return ids

    def _model_metabolite_ids(self, met):
        ids = {str(met.id)}
        stripped_id = self._strip_model_compartment(met.id, self._model_metabolite_compartment(met))
        if stripped_id:
            ids.add(stripped_id)

        for key in ["kegg.compound", "bigg.metabolite", "seed.compound"]:
            ids.update(self._annotation_values(met.annotation, key))

        return {value for value in ids if value}

    def _model_metabolite_output_id(self, met, m_name=None):
        if self.DB == "BIGG":
            if self.metabolite_id_compartments:
                return met.id
            return (
                self._select_annotation_value(met, "bigg.metabolite", m_name)
                or self._strip_model_compartment(met.id, self._model_metabolite_compartment(met))
                or met.id
            )
        if self.DB == "SEED":
            return (
                self._select_annotation_value(met, "seed.compound", m_name)
                or self._strip_model_compartment(met.id, self._model_metabolite_compartment(met))
                or met.id
            )
        if self.DB == "KEGG":
            return (
                self._select_annotation_value(met, "kegg.compound", m_name)
                or self._strip_model_compartment(met.id, self._model_metabolite_compartment(met))
                or met.id
            )
        return met.id

    def _database_metabolite_output_id(self, met):
        direct_id = self._database_annotation_value(met)
        if direct_id:
            return direct_id

        mapped_id = self._mapped_model_metabolite_database_id(met)
        if mapped_id:
            return mapped_id

        return (
            self._strip_model_compartment(met.id, self._model_metabolite_compartment(met))
            or met.id
        )

    def _database_annotation_value(self, met):
        key_by_db = {
            "BIGG": "bigg.metabolite",
            "SEED": "seed.compound",
            "KEGG": "kegg.compound",
        }
        key = key_by_db.get(self.DB)
        if key is None:
            return None
        return self._select_annotation_value(met, key)

    def _mapped_model_metabolite_database_id(self, met):
        if self.DB not in {"BIGG", "SEED", "KEGG"}:
            return None

        kegg_ids = self._annotation_values(met.annotation, "kegg.compound")
        for kegg_id in kegg_ids:
            mapped_id = self._mapped_kegg_metabolite_database_id(kegg_id)
            if mapped_id:
                return mapped_id

        source_ids = {
            "bigg": self._annotation_values(met.annotation, "bigg.metabolite"),
            "seed": self._annotation_values(met.annotation, "seed.compound"),
        }
        source_ids["bigg"].append(
            self._strip_model_compartment(met.id, self._model_metabolite_compartment(met))
        )

        for source_db, values in source_ids.items():
            for value in values:
                if not value:
                    continue
                for entry in self.m_mapper.reverse_lookup(source_db, value):
                    mapped_id = self._mapped_kegg_metabolite_database_id(entry.kegg)
                    if mapped_id:
                        return mapped_id

        return None

    def _mapped_kegg_metabolite_database_id(self, kegg_id):
        if self.DB == "KEGG":
            return kegg_id

        mapped = self.m_mapper.get(kegg_id)
        if not mapped:
            return None
        if self.DB == "BIGG":
            return mapped.bigg
        if self.DB == "SEED":
            return mapped.seed
        return None

    def _select_annotation_value(self, met, key, m_name=None):
        values = self._annotation_values(met.annotation, key)
        if not values:
            return None

        if m_name is not None:
            preferred = self._map_metabolite_ids(m_name)
            for value in values:
                if value in preferred:
                    return value

        return values[0]

    def _model_metabolite_compartment(self, met):
        compartment = getattr(met, "compartment", None)
        return str(compartment) if compartment else None

    def _normalize_compartment_filter(self, compartment):
        if compartment is None:
            return None

        compartment = str(compartment).strip()
        return compartment or None

    def _reaction_matches_compartment_filter(self, reaction):
        if self.compartment_filter is None:
            return True

        return self.compartment_filter in self._model_reaction_compartments(reaction)

    def _model_reaction_compartments(self, reaction):
        compartments = getattr(reaction, "compartments", None)
        if compartments:
            return {str(compartment) for compartment in compartments if compartment}

        found = set()
        for metabolite in getattr(reaction, "metabolites", {}):
            compartment = self._model_metabolite_compartment(metabolite)
            if compartment:
                found.add(compartment)
        return found

    def _model_metabolite_matches_compartment_filter(self, metabolite):
        if self.compartment_filter is None:
            return True

        return self._model_metabolite_compartment(metabolite) == self.compartment_filter

    def _strip_model_compartment(self, metabolite_id, compartment):
        if compartment and metabolite_id.endswith(f"_{compartment}"):
            return metabolite_id[: -(len(compartment) + 1)]
        return metabolite_id

    def _reaction_metabolite_entry(self, entry):
        reaction_metabolite = {
            "bigg_id": entry["bigg_id"],
            "coefficient": entry["coefficient"],
        }
        if entry.get("compartment"):
            reaction_metabolite["compartment"] = entry["compartment"]
        return reaction_metabolite
    
    def _subtract_not_in_model_reactions(self, global_idxs, all_nodes, anti_rs, r2indx_dict):

        for r_name in anti_rs:

            all_nodes[global_idxs["reactions"][r_name]] = None
            all_nodes[global_idxs["multimarkers"][r_name]["in"]] = None
            all_nodes[global_idxs["multimarkers"][r_name]["out"]] = None

            r2indx_dict[r_name] = None

        return all_nodes, r2indx_dict
    
    def _subtract_not_in_model_metabolites(
        self,
        global_idxs,
        nodes,
        anti_ms,
        reactions=None,
        reaction_index=None,
    ):
        referenced_nodes = set()
        if reactions is not None and reaction_index is not None:
            for r_name, reaction in reactions.items():
                if reaction_index.get(r_name) is None:
                    continue
                for segment in reaction.get("segments", {}).values():
                    referenced_nodes.add(segment["from_node_id"])
                    referenced_nodes.add(segment["to_node_id"])

        for m in anti_ms:
            node_idx = global_idxs["metabolites"][m]
            if node_idx in referenced_nodes:
                continue

            nodes[node_idx] = None

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

                entry = {
                    "bigg_id": self._secondary_metabolite_output_id(met),
                    "name": met.name,
                    "coefficient": coef,
                }
                compartment = self._model_metabolite_compartment(met)
                if (
                    (
                        self.use_model_metabolite_ids
                        or self.use_database_secondary_metabolite_ids
                    )
                    and self.metabolite_id_compartments
                    and compartment
                ):
                    entry["compartment"] = compartment

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

    def _secondary_metabolite_output_id(self, met):
        if self.use_model_metabolite_ids:
            return self._model_metabolite_output_id(met)
        if self.use_database_secondary_metabolite_ids:
            return self._database_metabolite_output_id(met)
        return met.id
    
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
                        compartment=entry.get("compartment"),
                    )
                    all_nodes[max_node_idx] = node
                    r_desc[r_name]["segments"][seg_counter] = self._prepare_edge_dict(
                        max_node_idx,
                        in_mm_idx,
                    )
                    r_desc[r_name]["metabolites"].append(
                        self._reaction_metabolite_entry(entry)
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
                        compartment=entry.get("compartment"),
                    )
                    all_nodes[max_node_idx] = node
                    r_desc[r_name]["segments"][seg_counter] = self._prepare_edge_dict(
                        out_mm_idx,
                        max_node_idx,
                    )
                    r_desc[r_name]["metabolites"].append(
                        self._reaction_metabolite_entry(entry)
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

        lateral_offset = self._calc_secondary_lateral_offset(index, total)
        pos = (
            anchor_pos
            + side * reaction_dir * self.secondary_metabolite_distance
            + perp * lateral_offset
        )

        return pos.tolist()

    def _calc_secondary_lateral_offset(self, index, total):
        lane = -(total - 1) / 2.0 + index
        if lane == 0:
            lane = 0.5

        return lane * self.secondary_metabolite_spacing


    def _generate_secondary_metabolite_dict(self, bigg_id, name, pos, compartment=None):

        node = {
            "node_type": "metabolite",
            "bigg_id": bigg_id,
            "name": name,
            "node_is_primary": False,
            "x": pos[0],
            "y": pos[1],
            "label_x": pos[0] + self.metabolite_label_shift[0],
            "label_y": pos[1] + self.metabolite_label_shift[1],
        }
        if compartment:
            node["compartment"] = compartment
        return node
    
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

        canvas["width"] = x - min_x + self.canvas_margin_x
        canvas["height"] = y - min_y + self.canvas_margin_y

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
    
