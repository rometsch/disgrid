########################################
# Fiel data class representing
# one specific field of one quantity
# at a given time

class Field:
    def __init__(self, grid, data, time):
        self.grid = grid
        self.data = data
        self.time = time

