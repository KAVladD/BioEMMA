import numpy as np

from utils.metanetx_mapper import MetaNetXMapper

class EscherMapper:

    def __init__(self, 
                 metabolites: dict, 
                 reactions: dict,
                 markers_dist: int = 10,
                 scaling_factor: float = 4,
                 metabolite_label_shift: list = [10, 10],
                 reaction_label_shift: list = [10, 10],
                 database: str = "BIGG",
                 
                 axis_epsilon: float = 2,):
        
        self.m_mapper = MetaNetXMapper("resources/metabolite_mapping.tsv", "first")
        self.r_mapper = MetaNetXMapper("resources/reaction_mapping.tsv", "first")

        self.metabolites = metabolites
        self.reactions = reactions

        self.map_main_metabolites = []
        
        self.markers_dist = markers_dist
        self.factor = scaling_factor
        self.metabolite_label_shift = metabolite_label_shift
        self.reaction_label_shift = reaction_label_shift

        # add to init
        self.h_margin = 100
        self.w_margin = 100
        self.mm_dist_part = 0.3
        self.use_const_mm_dist = False
        self.mm_dist_const = 300
        self.axis_epsilon = axis_epsilon   # порог "на одной оси" (пиксели)
        self.axis_offset = 20    # отступ от оси метаболитов
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

    def build_map(self, cobra_model):

        escher_map = []

        description = self._generate_description("test")
        escher_map.append(description)

        model = {}

        # prepare all metabolites descriptions
        m_desc, m2indx_dict = self._prepare_elements_descriptions(self.metabolites, self._generate_metabolite_dict)

        # prepare all reactions descriptions with subnodes
        r_desc, r2indx_dict = self._prepare_elements_descriptions(self.reactions, self._genarate_reaction_dict)
        r_nodes, r2node_dict = self._prepare_reactions_nodes(self.reactions)

        # prepare multimarkers between reactions and metabolites
        r_mm_nodes, r2mm_node_dict = self._prepare_reactions_multimarkers(self.reactions)

        # compose all nodes
        global_nodes_idxs = self._make_global_idxs(m2indx_dict, r2node_dict, r2mm_node_dict)
        all_nodes = self._compose_nodes(global_nodes_idxs, m_desc, r_nodes, r_mm_nodes)

        # update edges
        r_desc = self._add_edges_to_reactions_descriptions(self.reactions, r_desc, global_nodes_idxs)

        # extract and prepare data from model
        cobra_model_metabolites, cobra_model_reactions = self._parse_model(cobra_model)
        all_nodes, r2indx_dict = self._subtract_not_in_model_reactions(global_nodes_idxs, r_nodes, all_nodes, cobra_model_reactions, r2indx_dict)
        all_nodes = self._subtract_not_in_model_metabolites(global_nodes_idxs, all_nodes, m_desc, cobra_model_metabolites)

        model["nodes"] = {i:j for i,j in all_nodes.items() if j}
        
        model["reactions"] = {r2indx_dict[r]: r_desc[r] for r in r_desc.keys() if r2indx_dict[r]}

        model["text_labels"] = self.text_labels
        model["canvas"] = self.canvas
        
        model["nodes"], model["reactions"] = self._multiply_positions(model["nodes"], model["reactions"])
        model["canvas"] = self._tune_canvas(model["nodes"], model["canvas"] )
        model["nodes"], model["reactions"] = self._align_nodes(model["nodes"], model["reactions"], model["canvas"])

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
        m_dict["label_y"] = str(float(pos[1]) + self.metabolite_label_shift[0])

        return m_dict
    
    def _genarate_reaction_dict(self, ids, reaction, name=None):

        if self.DB == "BIGG":
            id = ids["BIGG"]
        elif self.DB == "SEED":
            id = ids["SEED"]
        elif self.DB == "KEGG":
            id = ids["KEGG"]

        reaction_dict = {
            "name": name,
            "bigg_id": id,
            "reversibility": reaction["reversibility"],
            "gene_reaction_rule": "",
            "genes": [],
            "metabolites": [],
            "segments": {},
        }

        reaction_dict["label_x"] = str(float(reaction["position"][0]) + self.reaction_label_shift[0])
        reaction_dict["label_y"] = str(float(reaction["position"][1]) + self.reaction_label_shift[1])

        reaction_dict["metabolites"].extend([{"kegg_id": self.metabolites[m]["ids"]["KEGG"],
                                              "bigg_id": self.metabolites[m]["ids"]["BIGG"],
                                              "seed_id": self.metabolites[m]["ids"]["SEED"], 
                                              "coef": -1} 
                                             for m in reaction["substrates"].get("main", [])])
        
        reaction_dict["metabolites"].extend([{"kegg_id": self.metabolites[m]["ids"]["KEGG"],
                                              "bigg_id": self.metabolites[m]["ids"]["BIGG"],
                                              "seed_id": self.metabolites[m]["ids"]["SEED"],
                                              "coef": 1} 
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
            print(center)

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
        """Вычисляет позиции обоих мультимаркеров для реакции"""
        
        reaction_pos = np.array(reaction_data["position"], dtype=np.float64)
    
        substrates = reaction_data.get("substrates", {})
        products = reaction_data.get("products", {})
        
        # Собираем ВСЕ метаболиты для проверки общей оси
        all_mets = substrates.get("main", []) + products.get("main", [])
        all_positions = [np.array(self.metabolites[m]["position"], dtype=np.float64) for m in all_mets]
        
        # Проверяем общую ось
        common_axis_type, common_axis_value = self._check_metabolites_on_same_axis(all_positions)
        
        # Теперь считаем позиции с учётом общей оси
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
        Вычисляет позицию одного мультимаркера.
        Приоритет:
        1. Метаболиты на одной оси между собой → mm на этой оси, x/y от реакции
        2. Есть метаболит на оси с реакцией → mm на этой оси, перпендикуляр от opposite
        3. Fallback → центр масс
        """
        
        mets = metabolites.get("main", None)
        if not mets:
            return None
        
        positions = [np.array(self.metabolites[m]["position"], dtype=np.float64) for m in mets]
        
        # === 1. Если есть общая ось для всех метаболитов реакции ===
        if common_axis_type == "vertical":
            mm_x = reaction_pos[0]
            # Сдвиг в сторону СВОИХ метаболитов
            my_center_y = np.mean([p[1] for p in positions])
            offset_dir = np.sign(my_center_y - reaction_pos[1])
            offset_dir = offset_dir if offset_dir != 0 else 1
            mm_y = reaction_pos[1] + self.axis_offset * offset_dir
            return [mm_x, mm_y]
        
        if common_axis_type == "horizontal":
            mm_y = reaction_pos[1]
            # Сдвиг в сторону СВОИХ метаболитов
            my_center_x = np.mean([p[0] for p in positions])
            offset_dir = np.sign(my_center_x - reaction_pos[0])
            offset_dir = offset_dir if offset_dir != 0 else 1
            mm_x = reaction_pos[0] + self.axis_offset * offset_dir
            return [mm_x, mm_y]
        
        # === 2. Проверяем: есть ли метаболит на оси с реакцией? ===
        aligned_type, aligned_pos = self._find_aligned_metabolite(positions, reaction_pos)
        
        if aligned_type is not None:
            opposite_mets = opposite_metabolites.get("main", [])
            
            if aligned_type == "horizontal":
                # Метаболит на горизонтали с реакцией
                mm_y = reaction_pos[1]
                
                # Сдвиг в сторону СВОЕГО метаболита, не opposite
                offset_dir = np.sign(aligned_pos[0] - reaction_pos[0])
                offset_dir = offset_dir if offset_dir != 0 else 1
                mm_x = reaction_pos[0] + self.axis_offset * offset_dir
                
                return [mm_x, mm_y]

            if aligned_type == "vertical":
                # Метаболит на вертикали с реакцией
                mm_x = reaction_pos[0]
                
                offset_dir = np.sign(aligned_pos[1] - reaction_pos[1])
                offset_dir = offset_dir if offset_dir != 0 else 1
                mm_y = reaction_pos[1] + self.axis_offset * offset_dir
                
                return [mm_x, mm_y]
        
        # === 3. Fallback: центр масс ===
        return self._calc_multimarker_position_by_mass_center(reaction_pos, positions)
    
    def _check_metabolites_on_same_axis(self, positions):
        """
        Проверяет, лежат ли все метаболиты на одной оси.
        Возвращает (axis_type, axis_value) или (None, None)
        """
        if len(positions) < 2:
            return None, None
        
        ys = [p[1] for p in positions]
        xs = [p[0] for p in positions]
        
        # Проверка горизонтали
        if max(ys) - min(ys) < self.axis_epsilon:
            return "horizontal", np.mean(ys)
        
        # Проверка вертикали
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
                # Оба условия — выбираем по МЕНЬШЕМУ отклонению
                # Меньший dy = ближе к горизонтали, меньший dx = ближе к вертикали
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
        Находит координату (x или y) ближайшего opposite метаболита.
        coord_idx: 0 для x, 1 для y
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
        """Fallback: позиция через центр масс"""
        
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

    def _parse_model(self, model):

        rs = self._extract_model_reactions(model)
        ms = self._extract_model_metabolites(model)

        return ms, rs
    
    def _extract_model_reactions(self, model):

        rs = {"KEGG": set(), "BIGG": set(), "SEED": set()}

        for _, r in enumerate(model.reactions):

            for k in rs.keys():

                names = r.annotation.get(f"{k.lower()}.reaction")
                names = names if type(names) is list else [names]
            
                rs[k] |= set(names)

        return rs

    def _extract_model_metabolites(self, model):

        ms = {"KEGG": set(), "BIGG": set(), "SEED": set()}

        for _, r in enumerate(model.metabolites):

            for k in ms.keys():

                names = r.annotation.get(f'{k.lower()}.{"metabolite" if k=="BIGG" else "compound"}')
                names = names if type(names) is list else [names]
            
                ms[k] |= set(names)

        return ms
    
    def _subtract_not_in_model_reactions(self, global_idxs, r_nodes, nodes, model_reactions, r2indx_dict):

        for r, _ in r_nodes.items():

            keggs = set([r])
            biggs = set(self.r_mapper[r].bigg_all) if self.r_mapper.get(r) else set()
            seeds = set(self.r_mapper[r].seed_all) if self.r_mapper.get(r) else set()

            if not (keggs & model_reactions["KEGG"] or biggs & model_reactions["BIGG"] or seeds & model_reactions["SEED"]):
                nodes[global_idxs["reactions"][r]] = None

                nodes[global_idxs["multimarkers"][r]["in"]] = None
                nodes[global_idxs["multimarkers"][r]["out"]] = None

                r2indx_dict[r] = None

        return nodes, r2indx_dict
    
    def _subtract_not_in_model_metabolites(self, global_idxs, nodes, m_nodes, model_ms):

        for m in m_nodes.keys():

            keggs = set([m])
            biggs = set(self.m_mapper[m].bigg_all) if self.r_mapper.get(m) else set()
            seeds = set(self.m_mapper[m].seed_all) if self.r_mapper.get(m) else set()

            if not (keggs & model_ms["KEGG"] or biggs & model_ms["BIGG"] or seeds & model_ms["SEED"]):
                nodes[global_idxs["metabolites"][m]] = None

        return nodes

    # scaling and canvas 

    def _tune_canvas(self, nodes, canvas):

        x, y = 0, 0

        for i, compound_dict in nodes.items():

            if not compound_dict:
                continue

            if compound_dict["x"]:

                if x < float(compound_dict["x"]):
                    x = float(compound_dict["x"])
                if y < float(compound_dict["y"]):
                    y = float(compound_dict["y"])

        canvas["width"] = x + self.h_margin
        canvas["height"] = y + self.w_margin

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

        x, y = 0, 0
        n_x, n_y = 0, 0

        for node in nodes.values():

            if not node:
                continue

            x += float(node["x"])
            y += float(node["y"])
            n_x += 1
            n_y += 1

            if "label_x" in node.keys():

                x += float(node["label_x"])
                y += float(node["label_y"])
                n_x += 1
                n_y += 1

        for reaction in reactions.values():

            if not reaction:
                continue

            x += float(reaction["label_x"])
            y += float(reaction["label_y"])
            n_x += 1
            n_y += 1

        return x/n_x, y/n_y
    