class EscherMerger():

    def __init__(self, maps=None):

        self.reset()

        self.maps = list(maps) if maps is not None else []
        self.bbox_size = [1, 1]

    def reset(self):

        self.maps = []
        self.bbox_size = [1, 1]

    def append(self, emap):

        self.maps.append(emap)

    def merge_maps(self):

        assert (self.maps), "okey, what r u planing 2 merg?"

        self._update_bbox_size()

        n = len(self.maps) ** (1/2)
        n = int(n) if int(n) == n else int(n) + 1

        shifts = [[[(n/2 - n + j) + 0.5, (n/2 - n + i) + 0.5] for j in range(n)] for i in range(n)]

        for i, emap in enumerate(self.maps):
            self.maps[i] = self._shift_positions(emap, shifts[i//n][i%n])

        merged_map = self.maps[0]
        for i in range(1, len(self.maps)):
            merged_map = self._merge_two_maps(merged_map, self.maps[i])

        # merged_map[1]["canvas"]["width"] = n * self.bbox_size[0] * 1.1
        # merged_map[1]["canvas"]["height"] = n * self.bbox_size[0] * 1.1

        merged_map = self._tune_canvas(merged_map)

        return merged_map

    def _update_bbox_size(self):

        for emap in self.maps:

            emap_size = self._get_map_size(emap)

            if self.bbox_size[0] < emap_size[0]:
                self.bbox_size[0] = emap_size[0]
            if self.bbox_size[1] < emap_size[1]:
                self.bbox_size[1] = emap_size[1]

    def _get_map_size(self, emap):

        x, y = 0, 0
        nodes = emap[1]["nodes"]

        for _i, node in nodes.items():

            if node["x"]:

                if x < float(node["x"]):
                    x = float(node["x"])
                if y < float(node["y"]):
                    y = float(node["y"])

        return [x, y]
    
    def _shift_positions(self, emap, shift):

        nodes = emap[1]["nodes"]
        reactions = emap[1]["reactions"]

        dx = shift[0] * self.bbox_size[0]
        dy = shift[1] * self.bbox_size[1]

        for node in nodes.values():

            node["x"] = float(node["x"]) + dx
            node["y"] = float(node["y"]) + dy

            if "label_x" in node.keys():

                node["label_x"] = float(node["label_x"]) + dx
                node["label_y"] = float(node["label_y"]) + dy

        for reaction in reactions.values():

            reaction["label_x"] = float(reaction["label_x"]) + dx
            reaction["label_y"] = float(reaction["label_y"]) + dy

        return emap
    
    def _merge_two_maps(self, base_emap, emap):

        base_reactions = base_emap[1]["reactions"]
        base_nodes = base_emap[1]["nodes"]

        next_reaction_id = self._next_numeric_key(base_reactions)
        next_node_id = self._next_numeric_key(base_nodes)
        next_segment_id = self._next_segment_id(base_reactions)

        node_id_map = {}
        for key, node in emap[1]["nodes"].items():
            new_key = str(next_node_id)
            next_node_id += 1
            node_id_map[str(key)] = new_key
            base_nodes[new_key] = node

        for _key, reaction in emap[1]["reactions"].items():

            segments_buffer = {}
            for _s_key, segment in reaction["segments"].items():
                segment["from_node_id"] = node_id_map[str(segment["from_node_id"])]
                segment["to_node_id"] = node_id_map[str(segment["to_node_id"])]

                segments_buffer[str(next_segment_id)] = segment
                next_segment_id += 1

            reaction["segments"] = segments_buffer
            base_reactions[str(next_reaction_id)] = reaction
            next_reaction_id += 1

        return base_emap

    def _next_numeric_key(self, values):
        if not values:
            return 0
        return max(int(key) for key in values.keys()) + 1

    def _next_segment_id(self, reactions):
        segment_ids = [
            int(segment_id)
            for reaction in reactions.values()
            for segment_id in reaction.get("segments", {}).keys()
        ]
        if not segment_ids:
            return 0
        return max(segment_ids) + 1
    
    def _tune_canvas(self, emap):

        x_min, x_max, y_min, y_max = 0, 0, 0, 0
        
        nodes = emap[1]["nodes"]

        for _i, node in nodes.items():

            if node["x"]:

                if x_min > float(node["x"]):
                    x_min = float(node["x"])
                if x_max < float(node["x"]):
                    x_max = float(node["x"])

                if y_min > float(node["y"]):
                    y_min = float(node["y"])
                if y_max < float(node["y"]):
                    y_max = float(node["y"])

        w = x_max - x_min
        h = y_max - y_min

        emap[1]["canvas"]["x"] = x_min - 0.05 * w
        emap[1]["canvas"]["y"] = y_min - 0.05 * h
        emap[1]["canvas"]["width"] = w * 1.1
        emap[1]["canvas"]["height"] = h * 1.1

        return emap

if __name__ == "__main__":

    m = EscherMerger()
    m.merge_maps()

