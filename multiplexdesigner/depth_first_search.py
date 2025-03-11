


#Depth-first search (DFS) is an algorithm used to traverse or search through a data structure,
#such as a graph or tree. The fundamental idea behind DFS is that it explores as far down a 
#branch of the graph or tree as possible before backtracking to explore alternative branches. 
#This approach contrasts with breadth-first search, which explores all nodes at the current depth 
#level before moving on to the next. You can read about breadth-first search HERE.
#
#Think of depth-first search as a way of exploring a maze: you choose a path and follow it until 
#you hit a dead end. Once you reach the end of that path, you backtrack to the last decision point 
#and try a different route. This deep exploration makes DFS an excellent choice when the goal is to 
#fully explore a structure, especially one with many potential paths.
# https://www.datacamp.com/tutorial/depth-first-search-in-python?dc_referrer=https%3A%2F%2Fwww.ecosia.org%2F

def depth_first_search_with_pruning(node, depth, alpha, beta, maximizing_player):
    # Base case: reached maximum depth or leaf node
    if depth == 0 or is_leaf_node(node):
        return evaluate_node(node)
    
    if maximizing_player:
        max_eval = float('-inf')
        for child in get_children(node):
            eval = depth_first_search_with_pruning(child, depth - 1, alpha, beta, False)
            max_eval = max(max_eval, eval)
            alpha = max(alpha, eval)
            
            # Pruning step: beta cut-off
            if beta <= alpha:
                break  # Cut off this branch
        return max_eval
    else:
        min_eval = float('inf')
        for child in get_children(node):
            eval = depth_first_search_with_pruning(child, depth - 1, alpha, beta, True)
            min_eval = min(min_eval, eval)
            beta = min(beta, eval)
            
            # Pruning step: alpha cut-off
            if beta <= alpha:
                break  # Cut off this branch
        return min_eval


