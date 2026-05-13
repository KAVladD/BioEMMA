class EscherMerger():

    def __init__(self, maps=[]):

        self.reset()

        self.maps = maps
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

        x,y = 0,0
        nodes = emap[1]["nodes"]

        for i, node in nodes.items():

            if node["x"]:

                if x < float(node["x"]): x = float(node["x"])
                if y < float(node["y"]): y = float(node["y"])

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

        base_reactions_n = len(base_emap[1]["reactions"])
        base_nodes_n = len(base_emap[1]["nodes"])
        base_last_edges = base_emap[1]["reactions"][str(base_reactions_n-1)]["segments"]

        base_edges_n = int(sorted(base_last_edges.keys())[-1])

        for key, reaction in emap[1]["reactions"].items():
            
            segments_buffer = {}
            for s_key, segment in reaction["segments"].items():
                segment["from_node_id"] = str(int(segment["from_node_id"]) + base_nodes_n)
                segment["to_node_id"] = str(int(segment["to_node_id"]) + base_nodes_n)

                segments_buffer[str(int(s_key) + base_edges_n)] = segment

            reaction["segments"] = segments_buffer
            base_emap[1]["reactions"][str(int(key) + base_reactions_n)] = reaction

        for key, node in emap[1]["nodes"].items():

            base_emap[1]["nodes"][str(int(key) + base_nodes_n)] = node

        return base_emap
    
    def _tune_canvas(self, emap):

        x_min, x_max, y_min, y_max = 0, 0, 0, 0
        
        nodes = emap[1]["nodes"]

        for i, node in nodes.items():

            if node["x"]:

                if x_min > float(node["x"]): x_min = float(node["x"])
                if x_max < float(node["x"]): x_max = float(node["x"])

                if y_min > float(node["y"]): y_min = float(node["y"])
                if y_max < float(node["y"]): y_max = float(node["y"])

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

