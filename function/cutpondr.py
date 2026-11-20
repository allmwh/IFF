import re
import pandas as pd
import requests
import urllib3


class CutPONDR():

    _ALGORITHM_FIELDS = {
        "VLXT": "VLXT",
        "XL1_XT": "XL1",
        "CAN_XT": "CAN",
        "VL3-BA": "VL3",
        "VSL2": "VSL2",
    }

    def __init__(self, show_progress_window=False):
        """
        automatically send sequences to pondr.com, and get order/disorder region information

        show_progress_window: kept for backward compatibility (unused in requests mode)
        """
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        self.session = requests.Session()

    def close(self):
        return self.session.close()

    def cut(self, sequence, protein_name="aa", algorithm="VLXT"):
        self.protein_name = protein_name
        self.sequence = sequence
        self.algorithm = algorithm

        if self.algorithm not in self._ALGORITHM_FIELDS:
            raise ValueError(
                f"Unsupported algorithm '{self.algorithm}'. "
                f"Choose from {list(self._ALGORITHM_FIELDS.keys())}"
            )

        base_url = "https://www.pondr.com/"
        field_name = self._ALGORITHM_FIELDS[self.algorithm]

        files = {
            field_name: (None, "on"),
            "ProteinName": (None, self.protein_name),
            "AccessionCode": (None, ""),
            "Sequence": (None, self.sequence),
            "wcwraw": (None, "on"),
            "submit_result": (None, "Submit Query"),
        }

        response = self.session.post(base_url, files=files, verify=False, timeout=60)
        response.raise_for_status()

        pre_blocks = re.findall(r"<pre>(.*?)</pre>", response.text, re.DOTALL | re.IGNORECASE)
        raw_block = None
        for block in reversed(pre_blocks):
            stripped = block.strip()
            if stripped.startswith("Num"):
                raw_block = stripped
                break

        if raw_block is None:
            raise RuntimeError("Failed to parse PONDR raw output from response.")

        self.result = raw_block
        lines = [line.strip() for line in raw_block.splitlines() if line.strip()]
        rows = [line.split() for line in lines[1:]]
        df = pd.DataFrame(rows, columns=["Num", "Res", "algorithm"])

        # make sequence_mask
        threshold = 0.5
        sequence_mask = ""
        for _, row in df.iterrows():
            num = float(row["algorithm"])

            if num >= threshold:
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
        same as get_sequence_mask(), but with different order/disorder symbols
        "0": order
        "1": disorder

        return: str
        """
        od_ident = self.__sequence_mask.replace("*", "0").replace(".", "1")
        return od_ident
