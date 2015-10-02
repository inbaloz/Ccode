import math

RADIUS_TO_STDEV = 0.75
DELTA = 0.01

class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Boundaries(object):
    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        

def grid_gaussian_overlap(center1, center2, radius1, radius2):
    
    # 1. Creating the gaussians.
    gaussian1 = create_gaussian_function(center1, radius1)
    gaussian2 = create_gaussian_function(center2, radius2)

    # 2. Calculate boundaries.
    boundaries = calculate_boundaries(center1, center2, radius1, radius2)

    # 3. Iterate over the grid and sum the overlap.
    return sum_grid_overlap(gaussian1, gaussian2, boundaries)

def create_gaussian_function(point, radius):
    stdev = RADIUS_TO_STDEV * radius
    def gaussian(x, y):
        exponent = - ( (x - point.x)**2 + (y - point.y)**2 ) / ( 2 * (stdev**2))
        return (1/(2*math.pi*(stdev**2))) * math.exp(exponent)

    return gaussian

def calculate_boundaries(point1, point2, radius1, radius2):
##    x_max = max(point1.x + radius1, point2.x + radius2)
##    x_min = min(point1.x - radius1, point2.x - radius2)
##    y_max = max(point1.y + radius1, point2.y + radius2)
##    y_min = min(point1.y - radius1, point2.y - radius2)
    
    x_max = max(point1.x + 10 * radius1, point2.x + 10 * radius2)
    x_min = min(point1.x - 10 * radius1, point2.x - 10 * radius2)
    y_max = max(point1.y + 10 * radius1, point2.y + 10 * radius2)
    y_min = min(point1.y - 10 * radius1, point2.y - 10 * radius2)

    return Boundaries(x_min, x_max, y_min, y_max)

def sum_grid_overlap(function1, function2, boundaries):
    overlap = 0

    curr_y = boundaries.y_min
    while curr_y <= boundaries.y_max:
        curr_x = boundaries.x_min
        while curr_x <= boundaries.x_max:
            overlap += (DELTA**2) * function1(curr_x, curr_y) * function2(curr_x, curr_y)
            curr_x += DELTA
        curr_y += DELTA

    return overlap


center1 = Point(0,0)
center2 = Point(0,0)
radius1 = 0.5
radius2 = 0.5
print grid_gaussian_overlap(center1, center2, radius1, radius2)
