import os

# EXTRACT "BÖLEN"

# path = "/cta/users/ardacetin/globalRepair/damageseq/Damage-seqProtocolFile/" For Damageseq
# path = "/cta/users/ardacetin/globalRepair/chipseq/ChIP-seqProtocolFile/" For ChIPseq
# path = "/cta/users/ardacetin/globalRepair/repair/XR-seqProtocolFile/" For XRseq

path = "/cta/users/ardacetin/globalRepair/chipseq/ChIP-seqProtocolFile/"
level_1_dir_list = os.listdir(path)

for i in level_1_dir_list:
    # only change "C:/Users/Versilov/Desktop/rpkm_study/chipseq/Chip-seqProtocolFile/"
    filelist = os.listdir(path + "{}".format(i))
    sorted_filelist = sorted(filelist)                           # yeni eklendi
    # change "6" to "?" for damageseq, "10" for repair
    denominator_file_name = sorted_filelist[9]                          # satır 15 e göre düzenlendi
    main_file_name = sorted_filelist[8]    
    # open the file with index "n"
    file_nth = open(path + "{}/{}".format(i, denominator_file_name))
    denominator_list = file_nth.readlines()
    denominator = int(denominator_list[0])
    file_nth_2 = open(path + "{}/{}".format(i, main_file_name))
    # devamı önceden yaptığımız kısım
    newlined_list_with_each_line = file_nth_2.readlines()
    list_with_each_line = []

    for g in newlined_list_with_each_line:
        list_with_each_line.append(g.replace("\n", ""))

    rpkm_list = []
    for h in list_with_each_line:
        split = h.split("\t")
        int_rpkm = ((int(split[6])) * (10**9))/(5*denominator)
        rpkm_list.append(str(int_rpkm))

    k = 0

    with open(path + "{}/".format(i) + "99_{}_RPKM_added.bed".format(i), 'w+') as file:
        for m in range(len(list_with_each_line)):
            file.write("{}\t{}\n".format(list_with_each_line[k],rpkm_list[k]))
            k += 1
