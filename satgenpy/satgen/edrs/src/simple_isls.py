
from itertools import combinations # gives out every pair of combination with no duplicates or (a,a)

def generate_leo_geo_isls(leo_ids, geo_ids):
    isls = []
    # LEO from index 0 to len(leo_ids)-1
    num_leo = len(leo_ids)
    
    # GEO from num_leo to num_leo + len(geo_ids)-1
    num_geo = len(geo_ids)
    
    for i in range(num_leo):
        for j in range(num_leo, num_leo + num_geo):
            isls.append((i, j))
            
    return isls

def generate_leo_leo_isls(leo_ids):
    # Idex 0-based for LEO
    leo_indices = list(range(len(leo_ids)))
    return list(combinations(leo_indices, 2))  # sorted before to be order from smaller to bigger so lists are (samll, big) as geo-leo

