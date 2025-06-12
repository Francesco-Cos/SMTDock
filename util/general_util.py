def generate_distance(p1, p2):
        diff = [e1 - e2 for e1, e2 in zip(p1, p2)]
        return sum([d**2 for d in diff])

def int_coord_distance(p1, p2, spacing):
        diff = [(e1 - e2)//spacing for e1, e2 in zip(p1, p2)]
        return diff

def center_of_mass(points):
    com = [0,0,0]
    for point in points:
        for i, coords in enumerate(point):
            com[i] += coords
    com = [el/len(points) for el in com]
    return com
