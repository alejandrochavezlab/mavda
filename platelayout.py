# -*- coding: utf-8 -*-

"""
--------------------------------------------------------------------------------
Copyright 2023 Benjamin Alexander Albert [Alejandro Chavez Lab at UCSD]
All Rights Reserved
MAVDA Academic License
platelayout.py
--------------------------------------------------------------------------------
"""

from typing import Tuple


class PlateLayout:

    @staticmethod
    def _getControlWells(plate : str) -> \
        Tuple[
            list[str],
            list[str],
            dict[str,list[str]]]:

        if (plate.startswith("NCIDIV_4861")   or \
            plate.startswith("NCIDIV_4862")   or \
            plate.startswith("NCIDIV_4863")   or \
            plate.startswith("NCIDIV_4864")   or \
            plate.startswith("NCIDIV_4869")   or \
            plate.startswith("NCIDIV_4870")   or \
            plate.startswith("NCIDIV_4871")   or \
            plate.startswith("NCIDIV_4872")   or \
            plate.startswith("HYCPK_1150" )) and \
           (('-'    not in plate) and \
            ("rep" not in plate)):
    
            negctrl = \
                ["D03","D04","D05","D06","D07",
                 "E03","E04","E05","E06","E07"]
            posctrl = ["D08","D09","D10","E08","E09","E10"]
            lopctrl = {"D10":["HIV-1","HIV-2"]}
        elif plate.startswith("NCIDIV_4865"):
            negctrl = \
                ["A07","B05","C08","D02","D09",
                 "E03","F11","G09","H03","H08"]
            posctrl = ["A03","B07","C04","E10","F04","G05"]
            lopctrl = {"A03":["HIV-1","HIV-2"]}
    
        elif plate.startswith("NCIDIV_4866"):
            negctrl = \
                ["A06","B10","C04","D03","D08",
                 "E02","F08","G09","H05","H11"]
            posctrl = ["A11","B05","C09","E07","F03","G04"]
            lopctrl = {"A11":["HIV-1","HIV-2"]}
    
        elif plate.startswith("NCIDIV_4867"):
            negctrl = \
                ["A08","B04","C10","D07","D09",
                 "E11","F04","G09","H03","H10"]
            posctrl = ["A05","B09","C04","E05","F06","G03"]
            lopctrl = {"A05":["HIV-1","HIV-2"]}
    
        elif plate.startswith("NCIDIV_4868") and '-' not in plate:
            negctrl = \
                ["A11","B09","C12","D01","D11",
                 "E10","F12","G09","H03","H10"]
            posctrl = ["A03","B03","C05","E03","F06","G04"]
            lopctrl = {"A03":["HIV-1","HIV-2"]}
    
        else:
            negctrl = [
                "A01","A12",
                "B01","B12",
                "C01","C12",
                "D12",
                "E01",
                "F01","F12",
                "G01","G12",
                "H01","H12"]
            posctrl = ["D01","E12"]
            lopctrl = {
                "D01":["HIV-1", "HIV-2"],
                "E12":["HIV-1", "HIV-2"]}
    
        return negctrl, posctrl, lopctrl

    @staticmethod
    def _getCompoundWellMap(plate : str) -> dict[str,str]:

        wells = ["{}{:02d}".format(let,num) for \
            let in list("ABCDEFGH") for num in range(1,13)]
        wellmap = {well:well for well in wells}

        def _swp(w1,w2):
            wellmap[w1] = w2
            wellmap[w2] = w1

        if plate.startswith("NCIDIV_4861") or \
           plate.startswith("NCIDIV_4862") or \
           plate.startswith("NCIDIV_4863") or \
           plate.startswith("NCIDIV_4864") or \
           plate.startswith("NCIDIV_4869") or \
           plate.startswith("NCIDIV_4870") or \
           plate.startswith("NCIDIV_4871") or \
           plate.startswith("NCIDIV_4872"):

            _swp("A01", "D03")
            _swp("B01", "D04")
            _swp("C01", "D05")
            _swp("D01", "D06")
            _swp("E01", "D07")
            _swp("F01", "D08")
            _swp("G01", "D09")
            _swp("H01", "D10")

            _swp("A12", "E03")
            _swp("B12", "E04")
            _swp("C12", "E05")
            _swp("D12", "E06")
            _swp("E12", "E07")
            _swp("F12", "E08")
            _swp("G12", "E09")
            _swp("H12", "E10")

        elif plate.startswith("NCIDIV_4865"):

            _swp("A01", "A03")
            _swp("B01", "B05")
            _swp("C01", "C04")
            _swp("D01", "D02")
            _swp("E01", "E03")
            _swp("F01", "F04")
            _swp("G01", "G05")
            _swp("H01", "H03")

            _swp("A12", "A07")
            _swp("B12", "B07")
            _swp("C12", "C08")
            _swp("D12", "D09")
            _swp("E12", "E10")
            _swp("F12", "F11")
            _swp("G12", "G09")
            _swp("H12", "H08")

        elif plate.startswith("NCIDIV_4866"):

            _swp("A01", "A06")
            _swp("B01", "B05")
            _swp("C01", "C04")
            _swp("D01", "D03")
            _swp("E01", "E02")
            _swp("F01", "F03")
            _swp("G01", "G04")
            _swp("H01", "H05")

            _swp("A12", "A11")
            _swp("B12", "B10")
            _swp("C12", "C09")
            _swp("D12", "D08")
            _swp("E12", "E07")
            _swp("F12", "F08")
            _swp("G12", "G09")
            _swp("H12", "H11")

        elif plate.startswith("NCIDIV_4867"):

            _swp("A01", "A05")
            _swp("B01", "B04")
            _swp("C01", "C04")
            _swp("D01", "D07")
            _swp("E01", "E05")
            _swp("F01", "F04")
            _swp("G01", "G03")
            _swp("H01", "H03")

            _swp("A12", "A08")
            _swp("B12", "B09")
            _swp("C12", "C10")
            _swp("D12", "D09")
            _swp("E12", "E11")
            _swp("F12", "F06")
            _swp("G12", "G09")
            _swp("H12", "H10")

        elif plate.startswith("NCIDIV_4868"):

            _swp("A01", "A03")
            _swp("B01", "B03")
            _swp("C01", "C05")
            # D01 is DMSO
            _swp("E01", "E03")
            _swp("F01", "F06")
            _swp("G01", "G04")
            _swp("H01", "H03")

            _swp("A12", "A11")
            _swp("B12", "B09")
            # C12 is DMSO
            _swp("D12", "D11")
            _swp("E12", "E10")
            # F12 is DMSO
            _swp("G12", "G09")
            _swp("H12", "H10")

        elif plate.startswith("HYCPK_1150"):

            _swp("A01", "D03")
            _swp("B01", "D04")
            _swp("C01", "D05")
            _swp("D01", "D06")
            _swp("E01", "D07")
            _swp("F01", "D08")
            _swp("G01", "D09")
            _swp("H01", "D10")

            _swp("A12", "E03")
            _swp("B12", "E04")
            _swp("C12", "E05")
            _swp("D12", "E06")
            _swp("E12", "E07")
            _swp("F12", "E08")
            _swp("G12", "E09")
            _swp("H12", "E10")

            _swp("H05", "F05")
            _swp("H06", "F06")
            _swp("H07", "F07")
            _swp("H08", "F08")
            _swp("H09", "F09")
            _swp("H10", "F10")
            _swp("H11", "F11")

        return wellmap


    def __init__(self, plate):
        self._plate = plate
        self._negctrl, self._posctrl, self._lopctrl = \
            PlateLayout._getControlWells(plate)
        self._compoundWellMap = \
            self._getCompoundWellMap(plate)
    
    def getPlate(self):
        return self._plate
    
    def getNegCtrl(self):
        return self._negctrl

    def getPosCtrl(self):
        return self._posctrl

    def getLopCtrl(self):
        return self._lopctrl

    def getCompoundWell(self, well : str) -> str:
        return self._compoundWellMap.get(well,well)
