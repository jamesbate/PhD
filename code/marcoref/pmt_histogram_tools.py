import data_analysis
from data_analysis import plot as plotta
import numpy as np
from data_analysis import fit
from scipy.special import factorial

def calculate_probability(n,t1,t2,tot):
		p0 = 0
		p1 = 0
		p2 = 0
		for i in range(len(n)):
			if i<t1:
				p0+=n[i]
			elif i>=t1 and i<t2:
				p1+=n[i]
			else:
				p2+=n[i]
		return p0/float(tot),p1/float(tot),p2/float(tot) #stupid python 2


def poissonian(x,mean):
	return (mean**x*np.exp(-mean))/(factorial(x))


is_sorted = lambda a: np.all(a[:-1] <= a[1:])

###TODO
# - Multiple scan DONE!
# - More than 2 ions
# - Upper limit histogram DONE!
# - Plot thresholds doesn't work DONE!
# - Histogram fit can be improved DONE!


class Pmt_histo_2ions(data_analysis.DefaultModule):
	#this function plots the histograms.
	#i just did it to understand the problem of the double click in the PMT
	#the problem did apear when we analized the parity flop data of 20121015
	module = "pmt_histo_analize"
	report_info = ""
	report_info_mini = ""

	
	def do(self, treshold1 =10 ,treshold2 = 45, plot_histogram = True,plot_parity = True,plot_treshold=True,fit_parity = True, frequency_parity = 1,fit_histogram = True):
		self.x = []
		self.y = []
		fit_results = []
		pmt_histogram = np.array([])

		histo = plotta.multipleplot(1)
		mp1 = plotta.multipleplot(1)
		parityplot = plotta.multipleplot(1)

		#fit functions
		cosine = fit.fit_function()
		cosine.name = "Cosine"
		cosine.detail = "AmplitudeP0 * np.cos(2*Pi*f*x + PhaseP1)"
		cosine.p = [None] * 2
		cosine.pdetail = ['Amplitude', 'Phase']
		cosine.pboundaries = [(0,1), (-np.pi,np.pi)]
		def getinits(x, y):
			return [abs(y.max()-y.mean()), 0]
		cosine.getinits = getinits
		def execute(x, *p):
			return p[0]*np.cos(2*np.pi*x*frequency_parity + p[1]) 
		cosine.execute = execute

		triple_poisson = fit.fit_function()
		triple_poisson.name = "Three Poissonians"
		triple_poisson.detail = "P0 * poisson(P1) + P2 * poisson(P3) + P4 * poisson(P5)"
		triple_poisson.p = [1,2,1,(treshold1+treshold2/2),1,2*treshold2]
		triple_poisson.pdetail = ['Amplitude1', 'Mean1','Amplitude2', 'Mean2','Amplitude3', 'Mean3']
		triple_poisson.pboundaries = [None,(0,treshold1),None,(treshold1,treshold2),None,(treshold2,None)]
		def getinits(x, y):
			return [1,2,1,(treshold1+treshold2/2),1,1.5*treshold2]
		triple_poisson.getinits = getinits
		def execute(x, *p):
			return p[0]*poissonian(x,p[1])+p[2]*poissonian(x,p[3])+p[4]*poissonian(x,p[5])
		triple_poisson.execute = execute




		plot2 = mp1.add_subplot()
		plot3 = parityplot.add_subplot()

		#For every scan selected
		for index_scan in range(self.number_of_scans):
				scan = self.scans[index_scan]
				cycles = int(scan.getCmdParameter('CYCLES'))
				qfp = scan.qfp
				print(qfp.file_fields)
				print(np.shape(qfp.file_fields))
				print(len(qfp.file_fields))


				scanned_variable = scan.getParameter("ExperimentCmd.scanSettings.dimension.1.ParamName")

				if scanned_variable !="":
					scanMode = scan.getParameter("ExperimentCmd.acquisitionMode")
					scanned_index = qfp.getHeaderIndex(scanned_variable)
					if scanned_index == -1:
						scanned_index = qfp.getHeaderIndex("transitions."+scanned_variable)
                
				print(scanned_variable)
				print(qfp.getColumn(scanned_index))
				scanned_variable = np.asarray(qfp.getColumn(scanned_index)).astype(np.float)
				# first_pmt_index = qfp.getHeaderIndex('addPMTcounts(1,0)')
				#scan_points = len(qfp.file_fields)
				print(scanned_variable)
				if is_sorted(scanned_variable):
					number_of_files = 1
				else:
					number_of_files = 2 
				scan_points = len(qfp.file_fields)/number_of_files
				# files = np.zeros(scan_points)
				#files = np.empty(shape=(number_of_files,scan_points,qfp.file_fields[first_pmt_index:-1].shape[1]))
				#print("Files: ",files)
				if number_of_files == 1:
					start = 0
					#first_pmt_index = qfp.getHeaderIndex('PMT1_1(0)')
				else:
					start = scan_points
				
				first_pmt_index = qfp.getHeaderIndex('PMT1_1(0)')
				
				x = qfp.file_fields[start:start+scan_points,first_pmt_index:]
				scanned_variable = qfp.getColumn(scanned_index)[start:start+scan_points]
				#	files[ii][:] = qfp.file_fields[scan_points*ii:scan_points*(ii+1)][first_pmt_index:-1]
				#	print files[ii]
				
				#first_pmt_index = qfp.getHeaderIndex('PMT1_2(0)')

				print number_of_files
				print first_pmt_index
				print cycles
				print scan_points
				#for i_number_of_files in range(scan_points-1):
					# for i_cycles in range(cycles):
						# index = int(i_number_of_files * cycles + i_cycles)
						# PMT_counts = qfp.file_fields[i_number_of_files+1][first_pmt_index + i_cycles]
						# pmt_histogram[i].append(PMT_counts)
				#x = qfp.file_fields
				length,cycles = np.shape(x)
				print(length,cycles)

				
				
				for j in range(length):
					pmt_histogram = np.concatenate((pmt_histogram,x[j]))
				

				upperbound_histo = 2*treshold2 + 1.5
				p0 = np.zeros(length)
				p1 = np.zeros(length)
				p2 = np.zeros(length)

				dp0 = np.zeros(length)
				dp1 = np.zeros(length)
				dp2 = np.zeros(length)

				for i in range(length):
					n,bins = np.histogram(x[i,:], np.arange(-0.5,upperbound_histo,1))
					p0[i],p1[i],p2[i] = calculate_probability(n,treshold1,treshold2,cycles)
					dp0[i] = np.max([np.sqrt(p0[i]*(1-p0[i])/cycles),1/float(cycles+2)])
					dp1[i] = np.max([np.sqrt(p1[i]*(1-p1[i])/cycles),1/float(cycles+2)])
					dp2[i] = np.max([np.sqrt(p2[i]*(1-p2[i])/cycles),1/float(cycles+2)])

				extralabel = ""	
				if self.number_of_scans>1:
					extralabel += " ["+scan.scan_number +"]"				
				plot2.errorbar(scanned_variable,p0,yerr = dp0,label="p0" + extralabel)
				plot2.errorbar(scanned_variable,p1,yerr = dp1,label="p1" + extralabel)
				plot2.errorbar(scanned_variable,p2,yerr = dp2,label="p2" + extralabel)
				plot2.legend()
				plot2.grid(True)

				
				parity = p0+p2-p1

				dparity = np.sqrt(dp0**2+dp1**2+dp2**2)

				extralabel = ""	
				if self.number_of_scans>1:
					extralabel += " ["+scan.scan_number +"]"

				plot3.errorbar(scanned_variable,parity,yerr = dparity,label = "parity"  + extralabel)
				plot3.grid(True)
				plot3.legend()

				if fit_parity:
					fit_func = cosine
					fit_func.p = [1,0]
					fit_results = fit.fit_utility(scanned_variable, parity, fit_func,yerr = dparity)
					print(fit_results.p)
					xgrid = np.linspace(np.min(scanned_variable),np.max(scanned_variable),num = 1000)
					plot3.plot(xgrid,fit_func.execute(xgrid,*fit_results.p),'r--')
					print(fit_func.detail)


				self.report_info+="<div id=\"myDIV"+str(index_scan) + "\" style=\"display: none;\">"
				if self.number_of_scans>1:
					self.report_info+="<h3>Scan"+" ["+scan.scan_number +"]"+ "</h3>"
				self.report_info+="<p>&lt;p0&gt;=" + str(np.average(p0,weights = 1/dp0**2)) + "+/-"+ str(np.sqrt(1/np.sum(1/dp0**2)))+ "</p>"
				self.report_info+="<p>&lt;p1&gt;=" + str(np.average(p1,weights = 1/dp1**2)) + "+/-"+ str(np.sqrt(1/np.sum(1/dp1**2)))+ "</p>"
				self.report_info+="<p>&lt;p2&gt;=" + str(np.average(p2,weights = 1/dp2**2)) + "+/-"+ str(np.sqrt(1/np.sum(1/dp2**2)))+ "</p>"
				self.report_info+="</div>"
				if fit_parity:
					self.report_info+="<div id=\"fit"+str(index_scan) + "\" style=\"display: none;\">"
					if self.number_of_scans>1:
						self.report_info+="<h3>Scan"+" ["+scan.scan_number +"]"+ "</h3>"
					self.report_info+="<p>Fit result: ("+ fit_func.detail + ")</br>" + str(fit_results.p) + "+/-"+ str(np.sqrt(np.diag(fit_results.covmat)))+ "</p>"
					self.report_info+="</div>"
				
				

		#PLT TOTAL HISTOGRAM
	    #np.arange doesn't include the last point
		plothisto = histo.add_subplot()	
		n, bins, patches = plothisto.hist(pmt_histogram, np.arange(-0.5,upperbound_histo,1), facecolor='g', alpha=0.75)
		if plot_treshold:
			limits = plothisto.get_ylim()
			plothisto.vlines(treshold1,*limits, color='r')
			plothisto.vlines(treshold2,*limits,color='r')

		bin_center = [ (bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]

		plothisto.grid(True)

		if fit_histogram:
			fit_func_poisson = triple_poisson
			fit_func_poisson.p = [100,2,100,(treshold1+treshold2/2),100,1.5*treshold2]
			fit_results_poisson = fit.fit_utility(bin_center,n,fit_func_poisson)
			xgrid = np.linspace(np.min(bin_center),np.max(bin_center),num = 1000)
			plothisto.plot(xgrid,fit_func_poisson.execute(xgrid,*fit_results_poisson.p),'b--')
			self.report_info+="<div id=\"fit_histo\" style=\"display: none;\">"
			self.report_info+="<p>Fit result: ("+ fit_func_poisson.detail + ")</br>" + str(fit_results_poisson.p) + "+/-"+ str(np.sqrt(np.diag(fit_results_poisson.covmat)))+ "</p>"
			self.report_info+="</div>"




		if plot_histogram:				
			self.plots.append(histo)
		self.plots.append(mp1)
		if plot_parity:
			self.plots.append(parityplot)

		self.report_info+="<button onclick=\"myFunction()\">Mean values</button><script>function myFunction(){var i; for (i = 0; i < "+str(self.number_of_scans)+ "; i++) {var x = document.getElementById(\"myDIV\" +i.toString());if (x.style.display === \"none\") {x.style.display = \"block\"; } else { x.style.display = \"none\";}}}</script>"
		self.report_info+="<button onclick=\"myFunction2()\">Fit results</button><script>function myFunction2(){var i; for (i = 0; i < "+str(self.number_of_scans)+ "; i++) {var x = document.getElementById(\"fit\"+i.toString()); var y = document.getElementById(\"fit_histo\"); if (x.style.display === \"none\") {x.style.display = \"block\"; } else { x.style.display = \"none\";} if (y.style.display === \"none\") {y.style.display = \"block\"; } else { y.style.display = \"none\";}}}</script>"


#class Pmt_histo_Nions(data_analysis.DefaultModule):
	#this function plots the histograms.
	#i just did it to understand the problem of the double click in the PMT
	#the problem did apear when we analized the parity flop data of 20121015
#	module = "pmt_histo_analize"
#	report_info = ""
#	report_info_mini = ""

	