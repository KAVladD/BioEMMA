import re

def make_m_id_dict(file):

    id_dict = {}

    kegg_regex = r'KEGG:(\s[CG]{1}\d{5};)*\s[CG]{1}\d{5}'
    kegg_compound = r'(\s[CG]{1}\d{5}){1}'

    bigg_regex = r'BiGG:( [\w\d]+;)* [\w\d]+'
    bigg_compound = r'([\w\d]+){1}'

    lines = file.readlines()

    for line in lines[1:-1]:#[1:-1]:

        f = False

        kegg = re.search(kegg_regex, line)
        bigg = re.search(bigg_regex, line)

        if not kegg:
            #counter[0] += 1
            f = True

        if not bigg:
            #counter[1] += 1
            f = True

        if f:
            continue

        kegg = kegg[0]
        kegg = re.findall(kegg_compound, kegg)
        kegg = [i.strip() for i in kegg]

        bigg = bigg[0][4:]
        bigg = re.findall(bigg_compound, bigg)
        bigg = [i.strip() for i in bigg]

        for bigg_id in bigg:
            for kegg_id in kegg:
                id_dict[bigg_id] = kegg_id

    return id_dict


def make_metabolites_id_dict():

    with open("resources/bigg_models_metabolites.txt") as f:

        d = make_bigg_m_id_dict(f)

    return d

def make_r_id_dict(file):

    id_dict = {}

    kegg_regex = r'KEGG:(\s[R]{1}\d{5};)*\s[R]{1}\d{5}'
    kegg_compound = r'(\s[R]{1}\d{5}){1}'

    bigg_regex = r'BiGG:( [\w\d]+;)* [\w\d]+'
    bigg_compound = r'([\w\d]+){1}'

    lines = file.readlines()

    counter = [0, 0]
    i=0

    for line in lines[1:-1]:#[1:-1]:
        i+=1

        f = False

        kegg = re.search(kegg_regex, line)
        bigg = re.search(bigg_regex, line)

        if 5370 <= i <= 5380:
            print(bigg)

        #print(bigg)

        if not kegg:
            counter[0] += 1
            f = True

        if not bigg:
            counter[1] += 1
            f = True

        if f:
            continue

        kegg = kegg[0]
        kegg = re.findall(kegg_compound, kegg)
        kegg = [i.strip() for i in kegg]

        bigg = bigg[0][5:]
        bigg = re.findall(bigg_compound, bigg)
        bigg = [i.strip() for i in bigg]

        for bigg_id in bigg:
            id_dict[bigg_id] = kegg

    # print(len(lines[1:-1]) - counter[0])
    # print(len(lines[1:-1]) - counter[1])

    return id_dict


def make_reactions_id_dict():

    with open("resources/bigg_models_reactions.txt") as f:

        #print(f.read())

        d = make_bigg_r_id_dict(f)

    return d


def make_bigg_r_id_dict(file):

    id_dict = {}

    kegg_regex = r'http\:\/\/identifiers\.org\/kegg\.reaction\/(R\d{5})'
    kegg_compound = r'(\s[R]{1}\d{5}){1}'

    ec_num_regex = r'http\:\/\/identifiers\.org\/ec-code\/(\d+\.\d+\.\d+\.\d+)'

    bigg_regex = r'BiGG:( [\w\d]+;)* [\w\d]+'
    bigg_compound = r'([\w\d]+){1}'

    #print(file.read())
    lines = file.readlines()

    counter = [0, 0]

    for line in lines[1:]:#[1:-1]:

        f = False

        line_data = line.split("\t")
        bigg = line_data[0]
        potential_kegg = line_data[4]
        old_biggs = line_data[5]

        kegg = re.findall(kegg_regex, potential_kegg)
        ec_num = re.findall(ec_num_regex, potential_kegg)

        #print(old_biggs)

        if not kegg:
            counter[0] += 1
            f = True

            if ec_num:
                #print(ec_num)
                counter[1] += 1


        if f:
            continue

        id_dict[bigg] = kegg


    #print(counter, len(lines))

    return id_dict

def make_bigg_m_id_dict(file):

    id_dict = {}

    kegg_regex = r'http\:\/\/identifiers\.org\/kegg\.compound\/(C\d{5})'
    kegg_compound = r'(\s[R]{1}\d{5}){1}'

    ec_num_regex = r'http\:\/\/identifiers\.org\/ec-code\/(\d+\.\d+\.\d+\.\d+)'

    bigg_regex = r'BiGG:( [\w\d]+;)* [\w\d]+'
    bigg_compound = r'([\w\d]+){1}'

    seed_regex = r'http\:\/\/identifiers\.org\/seed\.compound\/(cpd\d{5})'

    #print(file.read())
    lines = file.readlines()

    counter = [0, 0]

    for line in lines[1:]:#[1:-1]:

        f = False

        line_data = line.split("\t")
        if len(line_data) < 6:
            continue
        bigg = line_data[1]
        potential_kegg = line_data[4]
        old_biggs = line_data[5]

        kegg = re.findall(kegg_regex, potential_kegg)
        ec_num = re.findall(ec_num_regex, potential_kegg)
        seed = re.findall(seed_regex, potential_kegg)

        #print(seed, kegg)

        #print(old_biggs)

        if not kegg:
            counter[0] += 1
            f = True

            if ec_num:
                #print(ec_num)
                counter[1] += 1


        if f:
            continue

        id_dict[bigg] = kegg

        if seed:
            id_dict[seed[0]] = kegg


    #print(counter, len(lines))

    return id_dict


#print(len(make_metabolites_id_dict()))
#print(len(make_reactions_id_dict()))