"""This class simply saves all of the physical constants I need
"""
##----------------------------CLASS DEFINITION--------------------------------##

class PhysicalConstants:
    def __init__(self):
        self.h_bar = 1.054571817e-34
        self.m_e = 9.10938356e-31
        self.m_p = 1.6726219e-27
        self.m_calc = 40*self.m_p
        self.k_B = 1.38064852e-23
        self.h = 6.62607004e-34
