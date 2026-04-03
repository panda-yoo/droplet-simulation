
class Particle:
    def __init__(self, i, x, y, theta, radius) -> None:
        self.id = i
        # wrapped (for simulation)
        self.x = x
        self.y = y
        # unwrapped (for tracking)
        self.x_unwrapped = x
        self.y_unwrapped = y
        self.theta = theta
        self.radius = radius
        self.v = 2.0
        self.influence = 500.0

    def __str__(self):
        return f'theta is {self.theta}'