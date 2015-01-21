__author__ = 'michaelmclellan'


class Block:

    num_blocks = 0

    def __init__(self, sequence, query_id, target_id):

        self.unique_id = Block.num_blocks
        self.query_id = query_id
        self.target_id = target_id
        self.sequence = sequence

        Block.num_blocks += 1

    # def __del__(self):
    #     class_name = self.__class__.__name__
    #     print class_name, "baleted"


