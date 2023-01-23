# Script to be run WITHIN CASA!
day = "day4"
band="c"

filename = open("/home/cira/ATCA/bin/source_fluxesdict.json", "r+")
srcs = json.load(filename)

for key in list(srcs.keys()):
    if band == "c":
        ms = f"c3487_{day}_cx.ms"
        spw = "0"
    elif band == "x":
        ms = f"c3487_{day}_cx.ms"
        spw = "1"
    elif band == "l":
        ms = f"c3487_{day}_l.ms"
        spw = ""
    
    uvmodelfit(
        vis=ms,
        niter=10,
        field=key,
        selectdata=True,
        spw="0",
        outfile=f"{key}_{band}_{day}.cl"
    )
    tbl = cl.open(f"{key}_{band}_{day}.cl")
    fit = cl.getcomponent(0)
    flux = fit["flux"]["value"][0]
    srcs[key][f"{band}"][f"{day}"] = [flux]
    src_dict = {f"{key}": {f"{band}": {f"{day}": flux}}}
    # src_dict[f"{key}"][f"{band}"][f"{day}"] = [flux]


    with open(f"/home/cira/ATCA/bin/{key}_dict.json", "w") as f: 
        json.dump(src_dict,f)


# obsinfo = listobs(f"{ms}")
# info_keys = list(obsinfo.keys())
# for key in info_keys:
#     if key.startswith("scan"):
#         continue
#     else:
#         obsinfo.pop(key)

# for key in list(obsinfo.keys()):
#     fieldname = obsinfo[key]["0"]["FieldName"]
#     timestamp = obsinfo[key]["0"]["BeginTime"]
#     fitflux = srcs[fieldname][f"{band}"][f"{day}"]
#     scan = obsinfo[key]["0"]["scanId"]

#     uvmodelfit(
#         vis=ms,
#         niter=10,
#         field=key,
#         selectdata=True,
#         spw="0",
#         scan = scan,
#         sourcepar = [fitflux, 0, 0],
#         outfile=f"{fieldname}_{band}_{day}_{scan}.cl"
#     )
#     tbl = cl.open(f"{fieldname}_{band}_{day}_{scan}.cl")
#     fit = cl.getcomponent(0)
#     flux = fit["flux"]["value"][0]
#     srcs[key][f"{band}"][f"{timestamp}"] = [flux]


