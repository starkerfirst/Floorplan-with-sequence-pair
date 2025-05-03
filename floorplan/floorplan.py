################################################
#ELEC6910D: Floorplan
#Author: YANG Bohan
#Algorithm: Simulated Annealing with Sequence Pair Representation
################################################

import os
import random
import matplotlib.pyplot as plt
import math
import copy
import time
import tqdm

verbose = False
seed = 754
random.seed(seed)

def Kahn_topological_sort(self, relation_table, direction: str):
    # computing in-degree of each node
    num_modules = len(relation_table)
    in_degree = {node: 0 for node in range(num_modules)}
    coords = {i: 0 for i in range(self.num_modules)}
    for i in range(num_modules):
        for j in relation_table[i][direction]:
            in_degree[j] += 1

    # keep those nodes with in-degree 0 in the queue
    queue = [node for node in range(num_modules) if in_degree[node] == 0]
    sorted_nodes = []

    # from source to sink
    # O(V+E) much faster than evaluate_slow
    total_length = 0
    while queue:
        current_node = queue.pop(0)
        if direction == "right":
            length = self.module_dict[current_node].get_width()
        else:
            length = self.module_dict[current_node].get_height()
        sorted_nodes.append(current_node)
        for neighbor in relation_table[current_node][direction]:
            coords[neighbor] = max(coords[neighbor], coords[current_node] + length) # update coordinates of children nodes
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)
        total_length = max(total_length, coords[current_node] + length)

    return total_length, coords

def shortest_path_faster(sequence_pair, constraint_graph, source, destination, direction: str):
    assert(direction in ["right", "top"])
    # initialize the distance and predecessor arrays
    distance = {node: float('inf') for node in constraint_graph}
    predecessor = {node: None for node in constraint_graph}
    distance[source] = 0

    # initialize the queue
    queue = [source]

    while queue:
        current_node = queue.pop(0)
        for neighbor in constraint_graph[current_node][direction]:
            weight = sequence_pair.module_dict[current_node].get_width() if direction == "right" else sequence_pair.module_dict[current_node].get_height()
            if distance[current_node] + (-weight) < distance[neighbor]:
                distance[neighbor] = distance[current_node] + (-weight)
                predecessor[neighbor] = current_node
                if neighbor not in queue:
                    queue.append(neighbor)

    if destination == "all":
        return distance
    else:
        return distance[destination]

class Module:
    def __init__(self, height, width):
        self.height = height
        self.width = width
        self.rotation = False

    def rotate(self):
        self.rotation = not self.rotation

    def get_dim(self):
        if self.rotation:
            return (self.width, self.height)
        else:
            return (self.height, self.width)
        
    def get_height(self):
        if self.rotation:
            return self.width
        else:
            return self.height
        
    def get_width(self):
        if self.rotation:
            return self.height
        else:
            return self.width

class SequencePair:
    def __init__(self, max_chipwidth, numblocks, h_w_list):
        self.numblocks = numblocks
        self.max_chipwidth = max_chipwidth
        self.module_dict = {}
        self.num_modules = 0
        self.x_coords = None
        self.y_coords = None
        self.cur_width = None
        self.cur_height = None

        for h_w in h_w_list:
            self.module_dict[self.num_modules] = Module(h_w[0], h_w[1])
            self.num_modules += 1
        self.module_dict[-1] = Module(0, 0) # destination node
        self.module_dict[-2] = Module(0, 0) # source node
        assert(self.num_modules == numblocks)
        
        # init from fully vertical (must meet width constraint)
        valid = False
        self.positive_stepline_seq = list(range(self.num_modules)) # 0..N-1 init
        self.negative_stepline_seq = list(reversed(range(self.num_modules))) # 0..N-1 init
        while not valid:  
            random.shuffle(self.positive_stepline_seq[:self.num_modules//2]) 
            random.shuffle(self.positive_stepline_seq[self.num_modules//2:]) 
            random.shuffle(self.negative_stepline_seq[self.num_modules//2:]) 
            random.shuffle(self.negative_stepline_seq[:self.num_modules//2])
            if self.evaluate()[0] <= self.max_chipwidth:
                valid = True

    def random_select_i_j(self):
        if self.cur_width != None and self.cur_height != None:
            # get a random pair of modules, with higher probability for outer modules
            # calculate distance from center for each module
            center_x = self.cur_width / 2
            center_y = self.cur_height / 2
            distances = []
            for i in range(self.num_modules):
                distance = ((self.x_coords[i] - center_x) ** 2 + (self.y_coords[i] - center_y) ** 2) ** 0.5
                distances.append((i, distance))

            # sort modules by distance from center (descending)
            distances.sort(key=lambda x: x[1], reverse=True)

            # select first module with higher probability from outer modules
            if random.random() < 0.3:  
                i = distances[random.randint(0, self.numblocks // 8)][0]  
                j = i
                while i == j:
                    j = distances[random.randint(0, self.numblocks // 8)][0]  
                return i, j
            if random.random() < 0.5:  
                i = distances[random.randint(0, self.numblocks // 3)][0]  
                j = i
                while i == j:
                    j = distances[random.randint(0, self.numblocks // 3)][0]  
                return i, j 

        i = random.randint(0, self.num_modules - 1)
        j = i
        while i == j:
            j = random.randint(0, self.num_modules - 1)
        
        return i, j

    # swap a random pair of modules in the positive sequence
    def m1_op(self):
        i, j = self.random_select_i_j()
        # swap the modules
        self.positive_stepline_seq[i], self.positive_stepline_seq[j] = self.positive_stepline_seq[j], self.positive_stepline_seq[i]

    # swap a random pair of modules in both sequences
    def m2_op(self):
        i, j = self.random_select_i_j()  # get a random pair of modules
        # swap the modules
        self.positive_stepline_seq[i], self.positive_stepline_seq[j] = self.positive_stepline_seq[j], self.positive_stepline_seq[i]
        self.negative_stepline_seq[i], self.negative_stepline_seq[j] = self.negative_stepline_seq[j], self.negative_stepline_seq[i]

    # rotate a randomly selected module
    def m3_op(self):
        # get a random module
        i = random.randint(0, self.num_modules - 1)
        # rotate the module
        self.module_dict[i].rotate()

    # evaluate the solution and construct horizontal and vertical constraint graphs
    # return the width and height of the floorplan, and coordinates of the modules
    def evaluate(self):
        start_1 = time.perf_counter()
        # construct module relation table (by adjacency list)
        module_relation_table = [{"right":[], "top":[]} for i in range(self.num_modules)]

        # pre-computed index table
        pos_idx = {module: idx for idx, module in enumerate(self.positive_stepline_seq)}
        neg_idx = {module: idx for idx, module in enumerate(self.negative_stepline_seq)}

        for a in range(self.num_modules):
            pos_a = pos_idx[a]
            neg_a = neg_idx[a]
            for b in range(a+1, self.num_modules):
                pos_b = pos_idx[b]
                neg_b = neg_idx[b]

                # a is right of b if a is after b in both sequences
                if pos_a > pos_b and neg_a > neg_b:
                    module_relation_table[b]["right"].append(a)

                # a is left of b if a is before b in both sequences
                if pos_a < pos_b and neg_a < neg_b:
                    module_relation_table[a]["right"].append(b)

                # a is above b if a is before b in positive sequence and after b in negative sequence
                if pos_a < pos_b and neg_a > neg_b:
                    module_relation_table[b]["top"].append(a)

                # a is below b if a is after b in positive sequence and before b in negative sequence
                if pos_a > pos_b and neg_a < neg_b:
                    module_relation_table[a]["top"].append(b)

        end_1 = time.perf_counter()
        
        if verbose:
            print(f'Time taken to construct module relation table: {end_1 - start_1:.5f} seconds')
        
        ##########################################################
        # compute the width of the floorplan

        # In fact, it is not needed to do transistive reduction to compute the width and height of the floorplan
        # Topological sort can follow the constraint while calculating the dependency length
        # Kahn's algorithm for topological sort
        start_2 = time.perf_counter()
        chip_width, x_coords = Kahn_topological_sort(self, module_relation_table, "right")

        ##########################################################
        # compute the height of the floorplan     
        chip_height, y_coords = Kahn_topological_sort(self, module_relation_table, "top")
        end_2 = time.perf_counter()
        if verbose:
            print(f'Time taken to compute width and height: {end_2 - start_2:.5f} seconds')

        self.x_coords = x_coords
        self.y_coords = y_coords
        self.cur_height = chip_height
        self.cur_width = chip_width

        return chip_width, chip_height, x_coords, y_coords
    
    # implementation following course ppt
    # But it takes minutes to finish evaluating one solution
    def evaluate_slow(self):
        # construct module relation table (by adjacency list)
        module_relation_table = {i: {"left":[], "right":[], "top":[], "bottom":[]} for i in range(self.num_modules)}
        finished_pair = []

        # pre-computed index table
        pos_idx = {module: idx for idx, module in enumerate(self.positive_stepline_seq)}
        neg_idx = {module: idx for idx, module in enumerate(self.negative_stepline_seq)}

        for a in range(self.num_modules):
            pos_a = pos_idx[a]
            neg_a = neg_idx[a]
            for b in range(self.num_modules):
                if a == b:
                    continue
                if (a, b) in finished_pair or (b, a) in finished_pair:
                    continue
                finished_pair.append((a, b))

                pos_b = pos_idx[b]
                neg_b = neg_idx[b]

                # a is right of b if a is after b in both sequences
                if pos_a > pos_b and neg_a > neg_b:
                    module_relation_table[b]["right"].append(a)
                    module_relation_table[a]["left"].append(b)

                # a is left of b if a is before b in both sequences
                if pos_a < pos_b and neg_a < neg_b:
                    module_relation_table[b]["left"].append(a)
                    module_relation_table[a]["right"].append(b)

                # a is above b if a is before b in positive sequence and after b in negative sequence
                if pos_a < pos_b and neg_a > neg_b:
                    module_relation_table[b]["top"].append(a)
                    module_relation_table[a]["bottom"].append(b)

                # a is below b if a is after b in positive sequence and before b in negative sequence
                if pos_a > pos_b and neg_a < neg_b:
                    module_relation_table[b]["bottom"].append(a)
                    module_relation_table[a]["top"].append(b)
        if verbose:
            print("Module relation table:")
            for i in range(self.num_modules):
                print(f"Module {i}: {module_relation_table[i]}")
        
        ##########################################################
        # construct horizontal constraint graphs (by adjacency list)
        horizontal_constraint_graph = copy.deepcopy(module_relation_table)

        # Floyd-Warshall algorithm to calculate transitive closure
        distance_matrix = [[float('inf')] * self.num_modules for _ in range(self.num_modules)]
        for i in range(self.num_modules):
            distance_matrix[i][i] = 0
            for j in module_relation_table[i]["right"]:
                distance_matrix[i][j] = -self.module_dict[i].get_width()
        
        for k in range(self.num_modules):
            for i in range(self.num_modules):
                for j in range(self.num_modules):
                    if distance_matrix[i][k] != float('inf') and distance_matrix[k][j] != float('inf'):
                        distance_matrix[i][j] = min(distance_matrix[i][j], distance_matrix[i][k] + distance_matrix[k][j])
        
        # construct reduced graph
        for i in range(self.num_modules):
            for j in module_relation_table[i]["right"]:
                # check for other paths
                has_other_path = False
                # no need to update distance, since we only care whether it is connected
                for k in range(self.num_modules):
                    if k != i and k != j and distance_matrix[i][k] != float('inf') and distance_matrix[k][j] != float('inf'):
                        has_other_path = True
                        break
                
                if has_other_path:
                    horizontal_constraint_graph[i]["right"].remove(j)
                    horizontal_constraint_graph[j]["left"].remove(i)
                

        # append source and destination nodes
        horizontal_constraint_graph[-1] = {"left":[], "right":[]} # -1 is destination node
        horizontal_constraint_graph[-2] = {"left":[], "right":[]} # -2 is source node
        for i in range(self.num_modules):
            if len(horizontal_constraint_graph[i]["right"]) == 0:
                horizontal_constraint_graph[i]["right"] = [-1]
                horizontal_constraint_graph[-1]["left"].append(i)
            if len(horizontal_constraint_graph[i]["left"]) == 0:
                horizontal_constraint_graph[i]["left"] = [-2] 
                horizontal_constraint_graph[-2]["right"].append(i)
            
        if verbose:
            print("Horizontal constraint graph:")
            for i in range(self.num_modules):
                print(f"Module {i}: {horizontal_constraint_graph[i]}")

        # calculate the x-coordinates of the modules (lower-left corner)
        x_coords = {i: 0 for i in list(range(self.num_modules))+[-1]} # including destination node
        distance_list = shortest_path_faster(self, horizontal_constraint_graph, -2, "all", "right")
        for i in list(range(self.num_modules))+[-1]:
            x_coords[i] = -distance_list[i]
        chip_width = x_coords[-1]
        assert(max(x_coords.values()) == chip_width)
        if verbose:
            print("X-coordinates:")
            for i in range(self.num_modules):
                print(f"Module {i}: {x_coords[i]}")
            print(f"Chip width: {chip_width}")

        ##########################################################
        # construct vertical constraint graphs (by adjacency list)
        vertical_constraint_graph = copy.deepcopy(module_relation_table)
        
        # Floyd-Warshall algorithm to calculate transitive closure
        distance_matrix = [[float('inf')] * self.num_modules for _ in range(self.num_modules)]
        for i in range(self.num_modules):
            distance_matrix[i][i] = 0
            for j in module_relation_table[i]["top"]:
                distance_matrix[i][j] = -self.module_dict[i].get_height()
        
        for k in range(self.num_modules):
            for i in range(self.num_modules):
                for j in range(self.num_modules):
                    if distance_matrix[i][k] != float('inf') and distance_matrix[k][j] != float('inf'):
                        distance_matrix[i][j] = min(distance_matrix[i][j], distance_matrix[i][k] + distance_matrix[k][j])
        
        # construct reduced graph
        for i in range(self.num_modules):
            for j in module_relation_table[i]["top"]:
                # check for other paths
                has_other_path = False
                # no need to update distance, since we only care whether it is connected
                for k in range(self.num_modules):
                    if k != i and k != j and distance_matrix[i][k] != float('inf') and distance_matrix[k][j] != float('inf'):
                        has_other_path = True
                        break
            
                if has_other_path:
                    vertical_constraint_graph[i]["top"].remove(j)
                    vertical_constraint_graph[j]["bottom"].remove(i)

        # append source and destination nodes
        vertical_constraint_graph[-1] = {"top":[], "bottom":[]} # -1 is destination node
        vertical_constraint_graph[-2] = {"top":[], "bottom":[]} # -2 is source node
        for i in range(self.num_modules):
            if len(vertical_constraint_graph[i]["top"]) == 0:
                vertical_constraint_graph[i]["top"] = [-1]
                vertical_constraint_graph[-1]["bottom"].append(i)
            if len(vertical_constraint_graph[i]["bottom"]) == 0:
                vertical_constraint_graph[i]["bottom"] = [-2] 
                vertical_constraint_graph[-2]["top"].append(i)
            
        if verbose:
            print("Vertical constraint graph:")
            for i in range(self.num_modules):
                print(f"Module {i}: {vertical_constraint_graph[i]}")
            
        # calculate the y-coordinates of the modules (lower-left corner)
        y_coords = {i: 0 for i in list(range(self.num_modules))+[-1]} # including destination node
        distance_list = shortest_path_faster(self, vertical_constraint_graph, -2, "all", "top")
        for i in list(range(self.num_modules))+[-1]:
            y_coords[i] = -distance_list[i]
        chip_height = y_coords[-1]
        assert(max(y_coords.values()) == chip_height)
        if verbose:
            print("Y-coordinates:")
            for i in range(self.num_modules):
                print(f"Module {i}: {y_coords[i]}")
                print(f"Chip height: {chip_height}")

        return chip_width, chip_height, x_coords, y_coords

def collect_data(file):
    h_w_list = []
    file_path = os.path.join('2-floorplan-v0', file)
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        chipwidth = int(lines[0].split(':')[1].strip())
        numblocks = int(lines[1].split(':')[1].strip())

        for i in range(2, 2 + numblocks):
            parts = lines[i].strip().split(':')
            dims = parts[1].strip().split()
            h_w_list.append((int(dims[1]), int(dims[0])))
    
    assert(len(h_w_list) == numblocks)
    return chipwidth, numblocks, h_w_list

def verify_solution(sequence_pair):
    x_coords = sequence_pair.x_coords
    y_coords = sequence_pair.y_coords
    for i in range(sequence_pair.num_modules):
        for j in range(i+1, sequence_pair.num_modules):
            # check if any vertex of module i is inside module j
            if (x_coords[i] < x_coords[j] + sequence_pair.module_dict[j].get_width() and
                x_coords[i] + sequence_pair.module_dict[i].get_width() > x_coords[j] and
                y_coords[i] < y_coords[j] + sequence_pair.module_dict[j].get_height() and
                y_coords[i] + sequence_pair.module_dict[i].get_height() > y_coords[j]):
                print(f"Module {i} and Module {j} overlap!")
                return 
            # check if any vertex of module j is inside module i
            if (x_coords[j] < x_coords[i] + sequence_pair.module_dict[i].get_width() and
                x_coords[j] + sequence_pair.module_dict[j].get_width() > x_coords[i] and
                y_coords[j] < y_coords[i] + sequence_pair.module_dict[i].get_height() and
                y_coords[j] + sequence_pair.module_dict[j].get_height() > y_coords[i]):
                print(f"Module {j} and Module {i} overlap!")
                return 
    print("No overlap detected!")

def plot_solution(sequence_pair, file):
    fig, ax = plt.subplots()
    ax.set_title('Floorplan')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')

    max_width = 0
    max_height = 0
    for i in range(sequence_pair.num_modules):
        x = sequence_pair.x_coords[i]
        y = sequence_pair.y_coords[i]
        width = sequence_pair.module_dict[i].get_width()
        height = sequence_pair.module_dict[i].get_height()
        rect = plt.Rectangle((x, y), width, height, edgecolor='black', facecolor='grey')
        ax.add_patch(rect)
        # ax.text(x + width / 2, y + height / 2, str(i), ha='center', va='center')

        if x + width > max_width:
            max_width = x + width
        if y + height > max_height:
            max_height = y + height

    plt.xlim(0, max_width)
    plt.ylim(0, max_height)
    plt.gca().set_aspect('equal', adjustable='box')
    # plt.grid()
    plt.savefig(f'fig/floorplan_{file}.png')
    # plt.show()

def export_solution(sequence_pair, num_modules):
    x_coords = sequence_pair.x_coords
    y_coords = sequence_pair.y_coords
    if not os.path.exists(f"output/{num_modules}"):
        os.mkdir(f"output/{num_modules}")
    file_path = os.path.join(f"output/{num_modules}", f'floorplan.txt')
    with open(file_path, 'w') as f:
        # Index : x_coordinate y_coordinate rotation
        for i in range(sequence_pair.num_modules):
            f.write(f"{i} : {x_coords[i]} {y_coords[i]} {int(sequence_pair.module_dict[i].rotation)}\n")
    
    print(f"Solution exported to {file_path}")

def main():
    input_dir = '2-floorplan-v0'
    # RL_model_dir = 'RL_64_chiplet_my'
    file_list = os.listdir(input_dir)
    file_list = sorted([file for file in file_list if file.endswith('.txt')], key=lambda x: int(x.split('.')[0]))

    for file in file_list:
        start = time.time()
        print(f'Using file: {file}')
        max_chipwidth, numblocks, h_w_list = collect_data(file)
    
        sequence_pair = SequencePair(max_chipwidth, numblocks, h_w_list)
        # SA algorithm
        T_SCALE = 5
        total_rounds = 15000
        best_cost = float('inf')
        best_sol = [0, 0, {}, {}]

        for round in tqdm.trange(total_rounds):
            sequence_pair_trial = copy.deepcopy(sequence_pair)
            # generate a random operation
            op = random.randint(1, 3)
            if op == 1:
                sequence_pair_trial.m1_op()
            elif op == 2:
                sequence_pair_trial.m2_op()
            else:
                sequence_pair_trial.m3_op()

            # evaluate the solution
            new_chip_width, new_chip_height, x_coords, y_coords = sequence_pair_trial.evaluate()
            new_cost = new_chip_height if new_chip_width < max_chipwidth else float('inf')
        
            # Temperature function
            # T(x) = 0.07 * (1-x)/(1+8*x)
            x = round / total_rounds
            T = 0.07 * (1-x)/(1+8*x) * T_SCALE # 0.1 -> 0
                
            # Acceptance probability function
            if new_cost <= best_cost:
                sequence_pair = sequence_pair_trial
                best_cost = new_cost
                best_sol = [new_chip_width, new_chip_height, x_coords, y_coords]
            else:
                prob = math.exp(-((new_cost - best_cost)/best_cost)/T)
                if random.random() < prob:
                    sequence_pair = sequence_pair_trial
                    best_cost = new_cost
                    best_sol = [new_chip_width, new_chip_height, x_coords, y_coords]

        verify_solution(sequence_pair)

        print(f'Final cost: {best_cost}')
        print(f'Final chip width: {best_sol[0]}, chip height: {best_sol[1]}')
        if verbose:
            print(f'Final sequence pair: {sequence_pair.positive_stepline_seq}, {sequence_pair.negative_stepline_seq}')
            print(f'Final coordinates: {best_sol[2]}, {best_sol[3]}')

        # compute utilization
        ideal_area = 0
        for i in range(sequence_pair.num_modules):
            ideal_area += sequence_pair.module_dict[i].get_width() * sequence_pair.module_dict[i].get_height()
        real_area = best_sol[0] * best_sol[1]
        print(f'Utilization: {ideal_area / real_area:.2%}')
        
        export_solution(sequence_pair, numblocks)
        plot_solution(sequence_pair, file.split('.')[0])

        end = time.time()
        print(f'Time taken: {end - start:.2f} seconds')
        # breakpoint()

if __name__ == "__main__":
    main()