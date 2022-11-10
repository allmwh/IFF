import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import ui as UI
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.support import expected_conditions as EC


class CutPONDR():
    def __init__(self, show_progress_window=False):
        """
        automatically sends sequences to pondr.com, and get order/disrder info

        show_progress_window: show chrome during sequence identification by PONDR
        """
        options = webdriver.ChromeOptions()
        if not show_progress_window:
            options.add_argument("--headless")
            options.add_argument("--no-sandbox")
            options.add_argument("--disable-dev-shm-usage")
        self.driver = webdriver.Chrome(ChromeDriverManager().install(), chrome_options=options)

    def close(self):
        return self.driver.close()

    def click(self, box, check):
        if check == 1:
            if box.is_selected():
                pass
            else:
                box.click()
        elif check == 0:
            if box.is_selected():
                box.click()
            else:
                pass

    def cut(self, sequence, protein_name="aa", algorithm="VLXT"):
        self.protein_name = protein_name
        self.sequence = sequence
        self.algorithm = algorithm

        base_url = "http://www.pondr.com/"
        self.driver.get(base_url)

        # define bottom control
        # algorithm
        self.VLXT = self.driver.find_element(By.NAME,"VLXT")
        self.XL1_XT = self.driver.find_element(By.NAME,"XL1")
        self.CAN_XT = self.driver.find_element(By.NAME,"CAN")
        self.VL3_BA = self.driver.find_element(By.NAME,"VL3")
        self.VSL2 = self.driver.find_element(By.NAME,"VSL2")

        # output
        self.gra = self.driver.find_element(By.NAME,"graphic")
        self.stat = self.driver.find_element(By.NAME,"stats")
        self.report = self.driver.find_element(By.NAME,"seq")
        self.rawoutput = self.driver.find_element(By.NAME,"wcwraw")

        # switcher
        switcher = {
            "VLXT": lambda: self.click(self.VLXT, 1),
            "XL1_XT": lambda: self.click(self.XL1_XT, 1),
            "CAN_XT": lambda: self.click(self.CAN_XT, 1),
            "VL3-BA": lambda: self.click(self.VL3_BA, 1),
            "VSL2": lambda: self.click(self.VSL2, 1),
        }

        # click bottom
        self.click(self.VLXT, 0)
        self.click(self.VSL2, 0)
        self.click(self.gra, 0)
        self.click(self.stat, 0)
        self.click(self.report, 0)
        self.click(self.rawoutput, 1)
        switcher[self.algorithm]()

        # fill protein name
        name = self.driver.find_element(By.NAME,"ProteinName")
        name.clear()
        name.send_keys(self.protein_name)

        # fill sequence
        seqform = self.driver.find_element(By.NAME,"Sequence")
        seqform.clear()
        seqform.send_keys(self.sequence)

        # submit
        self.submit = self.driver.find_element(By.NAME,"submit_result")
        self.submit.click()

        # wait until loaded
        UI.WebDriverWait(self.driver, 2).until(
            EC.presence_of_element_located((By.XPATH, "/html/body/center/h2/img"))
        )

        # raw output
        # result analysis
        result = self.driver.find_element(By.XPATH,"/html/body/pre[6]")
        result = result.text
        self.result = result
        df = pd.DataFrame(result.split("\n"))
        df = pd.DataFrame(
            df[0].str.split(pat=" ", n=2).tolist(), columns=["Num", "Res", "algorithm"]
        )
        df = df.drop(0)

        # make sequence_mask
        thershold = 0.5
        sequence_mask = ""
        for index, row in df.iterrows():
            order = int(row[0])
            num = float(row[2])

            if num >= thershold:
                sequence_mask = sequence_mask + "."
            else:
                sequence_mask = sequence_mask + "*"
        self.__sequence_mask = sequence_mask

    def get_disorder_sequence(self):
        """
        get disorder sequence identified by PONDR, residues that identified as order are masked by "*"

        return: str, masked sequence
        """
        disorder_sequence = ""
        for x, y in zip(self.__sequence_mask, self.sequence):
            if x == ".":
                disorder_sequence = disorder_sequence + y
            elif x == "*":
                disorder_sequence = disorder_sequence + "*"
        self.__disorder_sequence = disorder_sequence
        return self.__disorder_sequence

    def get_order_sequence(self):
        """
        get order sequence identified by PONDR, position that identified as disorder are masked by "*"

        return: str, masked sequence
        """
        order_sequence = ""
        for x, y in zip(self.__sequence_mask, self.sequence):
            if x == ".":
                order_sequence = order_sequence + "*"
            elif x == "*":
                order_sequence = order_sequence + y
        self.__order_sequence = order_sequence
        return self.__order_sequence

    def get_sequence_mask(self):
        """
        residues are order or disorder,
        "*": order
        ".": disorder

        return: str
        """
        return self.__sequence_mask

    def get_od_ident(self):
        """
        same as get_sequence_mask, instead different symbols of order/disorder
        "0": order
        "1": disorder

        return: str
        """
        od_ident = self.__sequence_mask.replace("*", "0").replace(".", "1")
        return od_ident
