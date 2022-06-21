import os
import configparser
import pandas as pd
from .features import IndividualCosts, PairwiseCosts
from multiply.util.exceptions import NoPrimerNameException


class IndividualCostFactory:
    def __init__(self, ini_path, result_dir):
        """
        Create IndividualCosts from a configuration file stored
        at `ini_path`

        params
            ini_path: str
                Path to .ini file containing information
                about individual primer costs.

        """

        # Ensure points to valid file, then set
        if not os.path.isfile(ini_path):
            raise FileNotFoundError(
                f"No genome collection found at {ini_path}. Check file path is correct."
            )
        self._ini_path = ini_path

        # Read config object
        self._config = configparser.ConfigParser()
        self._config.read(ini_path)

        # Set results directory
        self.result_dir = result_dir

    @staticmethod
    def create_cost(cost_name, csv_path, column, weight):
        """
        Create an instance of `IndividualCosts`


        """

        # Check if data exists, else return
        if not os.path.exists(csv_path):
            print(f"No data found at {csv_path}.")
            print("Skipping -- will not be included in cost function.")
            return

        # Load data
        df = pd.read_csv(csv_path)

        # Sanity checks
        if not "primer_name" in df.columns:
            raise NoPrimerNameException(
                f"No column `primer_name` found in {csv_path}; costs must be assigned to a primers."
            )
        if not column in df.columns:
            print(f"Cost column {column} not found in {csv_path}.")
            print(f"Found columns: {', '.join(df.columns)}.")
            print("Skipping -- will not be included in cost function.")
            return

        # Convert target colum to pandas series
        primer_values = df[column]
        primer_values.index = df["primer_name"]

        return IndividualCosts(
            cost_name=cost_name, primer_values=primer_values, weight=weight
        )

    def get_individual_costs(self):
        """
        Get a list of individual costs

        """

        indv_costs = []
        for section in self._config.sections():

            # Create cost
            indv_cost = self.create_cost(
                cost_name=section,
                csv_path=f"{self.result_dir}/{self._config.get(section, 'file')}",
                column=self._config.get(section, "column"),
                weight=self._config.getfloat(section, "weight"),
            )

            # Store, if cost was successfully created
            if indv_cost is not None:
                indv_costs.append(indv_cost)

        return indv_costs


class PairwiseCostFactory:
    def __init__(self, ini_path, result_dir):
        """
        Create IndividualCosts from a configuration file stored
        at `ini_path`

        params
            ini_path: str
                Path to .ini file containing information
                about individual primer costs.

        """

        # Ensure points to valid file, then set
        if not os.path.isfile(ini_path):
            raise FileNotFoundError(
                f"No genome collection found at {ini_path}. Check file path is correct."
            )
        self._ini_path = ini_path

        # Read config object
        self._config = configparser.ConfigParser()
        self._config.read(ini_path)

        # Set results directory
        self.result_dir = result_dir

    @staticmethod
    def create_cost(cost_name, csv_path, weight):
        """
        Create an instance of `IndividualCosts`


        """

        # Check if data exists, else return
        if not os.path.exists(csv_path):
            print(f"No data found at {csv_path}.")
            print("Skipping -- will not be included in cost function.")
            return

        # Load data
        df = pd.read_csv(csv_path, index_col=0)  # importantly, there is an index here

        return PairwiseCosts(cost_name=cost_name, primer_values=df, weight=weight)

    def get_pairwise_costs(self):
        """
        Get a list of individual costs

        """

        pairwise_costs = []
        for section in self._config.sections():

            # Create cost
            pairwise_cost = self.create_cost(
                cost_name=section,
                csv_path=f"{self.result_dir}/{self._config.get(section, 'file')}",
                weight=self._config.getfloat(section, "weight"),
            )

            # Store, if cost was successfully created
            if pairwise_cost is not None:
                pairwise_costs.append(pairwise_cost)

        return pairwise_costs
