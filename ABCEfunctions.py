##########################################################################
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################


import sqlite3
import sys
import pandas as pd


def get_next_asset_id(db, first_asset_id):
    max_id = pd.read_sql("SELECT MAX(asset_id) FROM assets", db).iloc[0, 0]
    if max_id == None:
        next_id = first_asset_id
    else:
        next_id = max_id + 1

    return next_id



