#!/usr/bin/env python
from json import load, dumps
from os.path import exists, join
from sys import argv
from random import randint, choice


def print_state(state):
    # Note that %i will appear w/column name 'JOBID' in actual squeue output.
    # this is because %i shows the array-id if it's an array job and what we
    # consider the regular job-id if it's not an array job.
    print("JOBID,STATE")
    for job_id in state:
        if 'array_ids' in state[job_id]:
            # this is an array job
            for array_id in state[job_id]['array_ids']:
                if state[job_id]['array_ids'][array_id] <= 0:
                    end_state = state[job_id]['endgame'][array_id]
                else:
                    end_state = 'RUNNING'

                print(f"{array_id},{end_state}")
        else:
            # this is a non-array job
            if state[job_id]['countdown'] <= 0:
                end_state = state[job_id]['endgame']
            else:
                end_state = 'RUNNING'

            print(f"{job_id},{end_state}")


def generate_output(job_ids):
    results = {}

    for job_id in job_ids:
        is_successful = choice([True, False])
        is_array_job = choice([True, False])

        if is_array_job:
            result = {'job_id': job_id}
            result['array_ids'] = {}
            result['endgame'] = {}

            for i in range(0, randint(5, 15)):
                array_id = "%s_%d" % (job_id, i)
                result['array_ids'][array_id] = randint(3, 7)
                result['array_ids'][array_id] = randint(3, 7)
                if is_successful:
                    # all array jobs must be successful
                    result['endgame'][array_id] = "COMPLETED"
                else:
                    # some jobs may succeed but some may fail
                    result['endgame'][array_id] = choice(
                        ['COMPLETED', 'FAILED'])
            results[job_id] = result
        else:
            result = {'job_id': job_id}
            result['countdown'] = randint(3, 7)
            result['endgame'] = choice(['COMPLETED', 'FAILED'])
            results[job_id] = result

    return results


def save_state(state, file_path):
    with open(file_path, 'w') as f:
        print(dumps(state, indent=2), file=f)


def load_state(file_path):
    with open(file_path, 'r') as f:
        return load(f)


if __name__ == "__main__":
    # "squeue -t all -j " f"{','.join(job_ids)} " "-o '%i,%T'"
    job_ids = argv[4].split(',')

    state_file_path = join("sequence_processing_pipeline", "scripts",
                           "my_state.json")

    state = generate_output(job_ids)

    if exists(state_file_path):
        state = load_state(state_file_path)
    else:
        state = generate_output(job_ids)

    print_state(state)

    for job_id in state:
        if 'array_ids' in state[job_id]:
            # this is an array job.
            for array_id in state[job_id]['array_ids']:
                state[job_id]['array_ids'][array_id] -= 1
        else:
            # this is a standard job.
            state[job_id]['countdown'] -= 1

    save_state(state, state_file_path)
