from math import sqrt, sin, cos, tan, pi, radians

def nozzle_divergent_length(throat_radius: float, expansion_ratio: float, throat_angle: float, exit_angle:float):
    exit_radius = sqrt(expansion_ratio) * throat_radius
    xp = .382 * throat_radius * sin(throat_angle)
    yp = throat_radius + .382 * throat_radius * (1 - cos(throat_angle))
    ye = exit_radius
    cot_th = tan(pi/2 - throat_angle)
    cot_ex = tan(pi/2 - exit_angle)
    a = (cot_ex - cot_th) / (2*(ye - yp))
    b = cot_th - 2 * a * yp
    c = xp - a * yp**2 - b * yp
    xe = a * ye**2 + b*ye + c
    return xe

actual = nozzle_divergent_length(throat_radius=.1, expansion_ratio=50, throat_angle=radians(30), exit_angle=radians(5))
expected = 1.85

print(f'Expected and actual nozzle length are {expected:.2f} and {actual:.2f} meters respectively')