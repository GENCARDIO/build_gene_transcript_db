

                #print(result)
                # if field == "gene_name":
                #     if not result in gene_name:
                #         gene_name.append(result)
                #     else:
                #         same_gene_name += 1
                # if field == "ID":
                #     # there are genes with same identifier but with _PAR_Y
                #     # at the final of the ID e.g. ENSG00000229232.6_PAR_Y
                #     # (that we will differenciate between them as they are
                #     # pseudoautosomal region on the Y chromosome
                #     if "PAR_Y" not in result.split(".")[1]:
                #         # not taking into account gene version
                #         result = result.split(".")[0]
                #     else:
                #         result = str(result.split(".")[0]) + "_PAR_Y"

                # gene_dict[field] = result
            # if feature == "CDS":
            #     gene_name = next(filter(lambda item: item.startswith("gene="), info), None)
            #     isoform  = next(filter(lambda item: item.startswith("Parent="), info), None)
            #     isoform = isoform.replace("Parent=rna-", "")
            #     print(gene_name)
    # print(feature_set)



                # if gene_id != id:
                #     print(f"gene_id : {gene_id}\nid = {id}")


    # print(same_gene_name)        

#             difference = info_set - set(whole_gene_inf)
#             for dif in difference:
#                 difference_global.add(dif)
#                 # if inf.split("=")[0] not in info_set:
#                 #     not_always_present.add(inf.split("=")[0])
#                 # info_set.add(inf.split("=")[0])
        


#         if feature == "exon":
#             whole_exon_inf = list()

#             for inf in info:
#                 whole_exon_inf.append(inf.split("=")[0])
#                 info_exon_set.add(inf.split("=")[0])
#             difference_exon = info_exon_set - set(whole_exon_inf)
#             for dif in difference_exon:
#                 difference_exon_global.add(dif)
#                 # if inf.split("=")[0] not in info_set:
#                 #     not_always_present.add(inf.split("=")[0])
#                 # info_set.add(inf.split("=")[0])

#         if feature == "CDS":
#             whole_CDS_inf = list()

#             for inf in info:
#                 whole_CDS_inf.append(inf.split("=")[0])
#                 info_CDS_set.add(inf.split("=")[0])
#             CDS_difference = info_CDS_set - set(whole_CDS_inf)
#             for dif in CDS_difference:
#                 difference_CDS.add(dif)
        

# print(f"features = {features}")
# for info_s in info_CDS_set:
#     print(info_s)

# print(f"not always present: {difference_CDS}")
