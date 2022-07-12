import json
import matplotlib.pyplot as plt
import ctypes

if __name__ == '__main__':
    #simulation parameters/look compare with work of kraut
    simulation_data = {}
    m = 0
    #170h pro Teilung
    simulation_data["ReproduceConnected"] = {"b":1/170,"m":m, "alpha":1}

    #40h pro Teilung
    simulation_data["ReproduceFree"] = {"b":1/40,"m":m, "alpha":1}

    #200h pro Tod
    simulation_data["DieConnected"] = {"d_C":1/200,"c_CC":0.0001,"c_CF":0, "alpha_C":1, "alpha_F":1}

    #200h pro Tod
    simulation_data["DieFree"] = {"d_F":1/200,"c_FF":0,"c_FC":0,"alpha_F":1, "alpha_C":1}

    simulation_data["k_max"] = 1000000
    simulation_data["t_max"] = 2000
    simulation_data["K"] = 1
    simulation_data["step_size"] = 1

    # (10,10)
    simulation_data["start"] = (100,0,0,0)

    simulation_data["model"] = "connections" #connections/Kraut
    if simulation_data["model"] == "Kraut":
        #12h pro switch
        simulation_data["Switch_CF"] = {"s_CF":0, "alpha":1}

        #8h pro switch
        simulation_data["Switch_FC"] = {"s_FC":1/5000, "alpha":1, "alpha_C":0.8}


    elif simulation_data["model"] == "connections":
        simulation_data["addConnections"] = [{"p":1/8, "alpha":1},{"p":0.05, "alpha":1},{"p":0.05, "alpha":1}]
        simulation_data["removeConnections"] = [{"p":1/12, "alpha":1},{"p":0.05, "alpha":1},{"p":0.001, "alpha":1}]

        simulation_data["maxConnections"] = 3


    with open("modelParameters.json", 'w') as f:
        json.dump(simulation_data, f)

    lib = ctypes.cdll['./cmake-build-debug/libsimulationLibrary.so']
    lib['main']()

    #load simulation results
    f = open('cpp_results.json', )
    data = json.load(f)
    f.close()

    #transpose results for plotting
    numberVStime = []
    labels = []
    for i in range(len(data[0])):
        numberVStime.append(list(map(lambda x:x[i],data)))

    #norm number of cells to one
    totalCells = []
    for i in range(len(data)):
        totalCells.append(sum(data[i]))

    res = []
    for i in range(len(data)):
        res.append([])
        for k in range(len(data[i])):
            if totalCells[i]==0:
                res[i]=[0]*len(data[i])
            else:
                res[i].append(data[i][k]/totalCells[i])

    #transpose
    fractionsVStime = []
    for i in range(len(res[0])):
        fractionsVStime.append(list(map(lambda x:x[i],res)))

    fig, ax = plt.subplots()
    plts = []
    for i in range(len(data[0])):
        tmp, = ax.plot(numberVStime[i], label = str(i))
        plts.append(tmp)

    for i in range(len(data[0])):
        legend = ax.legend(handles=plts)


    plt.savefig("./pdfs/tmp2.pdf")
    plt.show()
    #plt.plot(nums[0])
    #plt.show()
    #plt.plot(nums[1])
    #plt.show()

    # Python program to create
    # a pdf file
    fig, ax = plt.subplots()
    plts = []
    for i in range(len(data[0])):
        tmp, = ax.plot(fractionsVStime[i], label = str(i))
        plts.append(tmp)

    for i in range(len(data[0])):
        legend = ax.legend(handles=plts)

    plt.savefig("./pdfs/tmp3.pdf")
    plt.show()

    #create pdf
    from fpdf import FPDF
    pdf = FPDF()

    # Add a page
    pdf.add_page()

    # set style and size of font
    # that you want in the pdf
    pdf.set_font("Arial", size = 13)
    def st(x):
        if type(x)==float:
            return "{:.2f}".format(x)
        else:
            return str(x)

    # create cell
    h = 8
    for key,value in simulation_data.items():
        if type(value) == dict and len(value)>3:
            pdf.cell(200, h, txt = st(key)+": ",
                     ln = 1, align = 'L')
            for key2,value2 in value.items():
                pdf.cell(200, h, txt = "          "+st(key2)+": "+st(value2),
                         ln = 1, align = 'L')
        elif type(value) == list and len(value)>2:
            pdf.cell(200, h, txt = st(key)+": ",
                     ln = 1, align = 'L')
            for item in value:
                pdf.cell(200, h, txt = "          "+st(item),
                         ln = 1, align = 'L')
        else:
            pdf.cell(200, h, txt = st(key)+": "+st(value),ln = 1, align = 'L')

    pdf.output("./pdfs/tmp1.pdf")
    from PyPDF2 import PdfFileMerger, PdfFileReader
    merger = PdfFileMerger()
    merger. append(PdfFileReader(open("./pdfs/tmp1.pdf", 'rb')))
    merger. append(PdfFileReader(open("./pdfs/tmp2.pdf", 'rb')))
    merger. append(PdfFileReader(open("./pdfs/tmp3.pdf", 'rb')))
    from datetime import datetime
    from datetime import date

    today = date.today()
    current_date = today.strftime("%d/%m/%Y")
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")

    #uncomment this to save pdf
    #merger.write("./pdfs/"+current_date+"_"+current_time+".pdf")