# download ISH images from Allen Brain Atlas for genes in the input file
# Xiaoyan, 2017

# import urllib
import urllib2
import os
# from joblib import Parallel, delayed
# from multiprocessing import Pool, Process
# import multiprocess
from allensdk.api.queries.image_download_api import ImageDownloadApi
from allensdk.api.queries.ontologies_api import OntologiesApi

# load allensdk classes
ida = ImageDownloadApi()
oa = OntologiesApi()


# # functions
# def downimage(currentsection):
#     global dirname
#     ida.download_section_image(currentsection, file_path=dirname + "/" + currentsection + ".jpg", downsample=5, quality=0.5)
#     return int(currentsection)
#
#
# def addwrap(currentsection):
#     p = Pool(processes=4)
#     newnum = p.apply_async(downimage, [currentsection])
#     print newnum.get()


# get the list of genes need to be targeted
genes = []
with open(r"E:\PROOOJECTS\12_Neuron_mapping\Genes\genes.txt", 'r') as f:
    genes = [line.strip('\n') for line in f]
genes = list(set(genes))

# # query all data set id of ISH experimenets and save as CSV file
# sectionset_url = "http://api.brain-map.org/api/v2/data/query.csv?criteria=" \
#                  "model::SectionDataSet," \
#                  "rma::criteria,[failed$eqfalse],products[abbreviation$eq'Mouse']," \
#                  "treatments[name$eq'ISH']," \
#                  "genes,plane_of_section," \
#                  "rma::options," \
#                  "[tabular$eq'plane_of_sections.name+as+plane','genes.acronym+as+gene','data_sets.id+as+section_data_set_id']," \
#                  "[order$eq'plane_of_sections.name,genes.acronym,data_sets.id']" \
#                  "&num_rows=1000000&start_row=0"
# downloadfile = urllib.URLopener()
# downloadfile.retrieve(sectionset_url, "ISHSectionSets.csv")

# search section id, gene by gene and download down-sampled section images
sections = []
for gene in genes:
    print gene
    dirname = "E:/PROOOJECTS/12_Neuron_mapping/Experiments/161229_V1ProbesRe/Allen/LowRes/" + gene
    try:
        os.mkdir(dirname)
        # get coronal data set for the selected gene
        section_url = "http://api.brain-map.org/api/v2/data/query.csv?criteria=" \
                      "model::SectionDataSet," \
                      "rma::criteria," \
                      "[failed$eqfalse],products[abbreviation$eq'Mouse']," \
                      "treatments[name$eq'ISH']," \
                      "genes[acronym$eq'" + gene + "']," + \
                      "plane_of_section[name$eq'coronal'],section_images," \
                      "rma::options," \
                      "[tabular$eq'plane_of_sections.name+as+plane','genes.acronym+as+gene'," \
                      "'data_sets.id+as+section_data_set_id'," \
                      "'sub_images.id+as+section_image_id']" \
                      "[order$eq'plane_of_sections.name,genes.acronym,data_sets.id,sub_images.id']" \
                      "&num_rows=1000000&start_row=0"
        downloadfile = urllib2.urlopen(section_url)
        downloaded = downloadfile.read()
        sectionlist = downloaded.split('\n')

        # get sagittal data set if coronal set is missing
        if len(sectionlist) == 2:
            section_url = "http://api.brain-map.org/api/v2/data/query.csv?criteria=" \
                          "model::SectionDataSet," \
                          "rma::criteria," \
                          "[failed$eqfalse],products[abbreviation$eq'Mouse']," \
                          "treatments[name$eq'ISH']," \
                          "genes[acronym$eq'" + gene + "']," \
                          "plane_of_section[name$eq'sagittal'],section_images," \
                          "rma::options," \
                          "[tabular$eq'plane_of_sections.name+as+plane','genes.acronym+as+gene'," \
                          "'data_sets.id+as+section_data_set_id'," \
                          "'sub_images.id+as+section_image_id']" \
                          "[order$eq'plane_of_sections.name,genes.acronym,data_sets.id,sub_images.id']" \
                          "&num_rows=1000000&start_row=0"
            downloadfile = urllib2.urlopen(section_url)
            downloaded = downloadfile.read()
            sectionlist = downloaded.split('\n')

        sectiondownload = []
        sectionset = 0
        if len(sectionlist) > 1:
            for section in sectionlist[1:]:
                section = section.split(',')
                if len(section) > 1 and (sectionset == 0 or sectionset == int(section[2])):
                    sections.append(section)
                    sectionset = int(section[2])
                    sectiondownload.append(section[3])
                else:
                    break

        for section in sectiondownload:
            ida.download_section_image(section, file_path=dirname + "/" + section + ".jpg", downsample=5,
                                       quality=50)
            # TODO: parallelization
            # pool = multiprocess.Pool(processes=4)
            # res = pool.map_async(downimage, sectiondownload)
            # res.get(timeout=1)
            # pool.close()
            # pool.join()
            # Parallel(n_jobs=6)(delayed(ida.download_section_image(section, file_path=dirname + "/" + section + ".jpg", downsample=5, quality=0.5)) for section in sectiondownload)
            #
            # Process(target=addwrap, args=sectiondownload).start()
    except (WindowsError, urllib2.URLError):
        pass        # skip already existing ones, escape when connection error happens


# write downloaded files
with open("E:/PROOOJECTS/12_Neuron_mapping/Genes/AllenISH/" + "sections.csv", 'w') as f:
    for line in sections:
        f.write("%s,%s,%s,%s\n" % tuple(line))

# download original images
for gene in genes:
    print gene
    dirname = "E:/PROOOJECTS/12_Neuron_mapping/Genes/AllenISH/HighResCoronal/" + gene
    files = os.listdir(dirname)
    for newsection in files:
        if not (newsection[:4]=='High'):
            if not os.path.isfile(dirname + '/High_' + newsection):     # skip already existing ones
                sectionname = newsection.split('.jpg')
                ida.download_section_image(sectionname[0], file_path=dirname + "/High_" + sectionname[0] + ".jpg", downsample=0, quality=100)
            os.remove(dirname + '/' + newsection)

