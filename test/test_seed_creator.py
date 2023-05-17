##########################################################################
## Copyright 2023 Argonne National Laboratory
#
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


import pytest
import seed_creator

abce_db_path = "this/is/a/test/abce_db.db"

def test_ask_user_permission_to_delete_yes():
    """
    Tests user permission. Where the answer is yes. 
    """

    seed_creator.input = lambda _ : 'y'

    output = seed_creator.ask_user_permission_to_delete(abce_db=abce_db_path)

    assert output == True

    return


def test_ask_user_permission_to_delete_no():
    """
    Tests user permission. Where the answer is no. 
    """

    seed_creator.input = lambda _ : 'n'

    output = seed_creator.ask_user_permission_to_delete(abce_db=abce_db_path)

    assert output == False

    return


@pytest.mark.skip("Unsure how to test with while loop")
def test_ask_user_permission_to_delete_bad_input():
    """
    Tests user permission. Where the answer is invalid. 
    """

    seed_creator.input = lambda _ : 'B'
    output = seed_creator.ask_user_permission_to_delete(abce_db=abce_db_path)

    return
