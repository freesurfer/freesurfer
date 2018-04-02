def intracranial_volume(structures, includeStructures):
    # sum over included structure volumes
    return sum(structure['vol'] for structure in structures if structure['name'] in includeStructures)
