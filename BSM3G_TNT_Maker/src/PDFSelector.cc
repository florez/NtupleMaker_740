// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/PDFSelector.h"



#include <fstream>
#include <cstdlib>

PDFSelector::PDFSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC, std::vector<std::string>& pdf_name):baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the PDFSelector Constructor --> calling SetBranches()." << std::endl;
  // Current pdfsets for the Pdf set mapping
  lheProducerToken=iCC.consumes< LHERunInfoProduct, edm::InRun >(edm::InputTag("externalLHEProducer"));
  lheProductToken=iCC.consumes< LHEEventProduct >(edm::InputTag("externalLHEProducer"));

  pdf_sets=&pdf_name;

  _is_data                      = iConfig.getParameter<bool>("is_data");
  SetBranches();
}

PDFSelector::~PDFSelector(){
  delete tree_;
}

void PDFSelector::Fill(const edm::Event& iEvent){
  Clear();


  /// Get the PDF weights for calculating uncertainties in NLO samples and store them in the userrecord

  edm::Handle<LHEEventProduct> lheInfoHandel;
  iEvent.getByToken(lheProductToken , lheInfoHandel);

  if (lheInfoHandel.isValid()) {//Begin of LHE PDF storing code
    ///We have to take care of the fact that in the weights vector, there are at first the scale variation weights and then the pdf variation weights.
    ///We put the scale weigths into a string (since the userrecord will contain them all as a string anyway)
    ///For pdf weigths, we deal with a vector because we have to sort the weights into the pdf groups using the nice map we created in beginRun()


    //sthe pdf weights vector has the same order as the pdf_sets and pdf_ids vector

    for (unsigned int i=0; i<lheInfoHandel->weights().size();i++){
      pdf_weights.push_back(lheInfoHandel->weights()[i].wgt);
    }

  }//End of LHE PDF storing code


}

void PDFSelector::SetBranches(){
  if(debug_) std::cout << "     PDFSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&pdf_weights               ,"pdf_weights");

  if(debug_) std::cout << "     PDFSelector: Finished setting branches." << std::endl;
}

void PDFSelector::Clear(){
  pdf_weights.clear();
}


void PDFSelector::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup){

  // In the begin of the run, we have to save all the pdfset relevant informations, in case our sample is NLO, so we can calculate the uncertainties later
  /// First, we create the list of pdf set identifyers which are used to determine the pdfset used during the production of the MC sample

  edm::Handle<LHERunInfoProduct> run;
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;


  //this will produce a warning but I found no other way to make it work!
  iRun.getByLabel("externalLHEProducer", run );
  //iRun.getByToken(lheProducerToken, run );

  if (run.isValid()) {//Begin of LHE PDF stuff preparing
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());
    ///Read them from the input file, people are welcome to improve this (make it work on the grid? :) )
    //ifstream infile(PdfSetFileName_.fullPath(),ios::in);
    //int number=0;  //Variable to hold each number as it is read

    int start = 0;
    int length = 0;


    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
      std::vector<std::string> lines = iter->lines();

      bool reading_scale = false;
      bool reading_pdfvar = false;

      string id_delimeter = "id=";
      string scale_muR_delimeter = "muR=";
      string scale_muF_delimeter = "muF=";
      string endofline_delimeter = "</weight>";
      string pdf_delimeter ="";
      string current_pdf="";
      string find_pdf="type=\"";


      string muR = "";
      string muF = "";

      for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
        string line=lines.at(iLine);
        //cout<<line;
        ///After getting the line, we set the corresponding bools to the let the skimmer know what to save
        if (line.find("envelope")!=string::npos && (line.find("Central scale variation")!=string::npos or line.find("scale_variation")!=string::npos)){
          if(line.find("Central scale variation")!=string::npos){
            scale_muR_delimeter = "mur=";
            scale_muF_delimeter = "muf=";
          }
          reading_scale = true;
          continue;
        }
        if (line.find("combine=\"hessian\"")!=string::npos  or line.find("PDF_variation")!=string::npos){
          start = line.find(find_pdf);
          if(line.find(find_pdf)==string::npos){
            find_pdf="name=\"";
          }
          length = line.find("\">")-start ;
          current_pdf=line.substr(start+find_pdf.size(), length-find_pdf.size());
          reading_pdfvar = true;
          continue;
        }
        ///This is the way every group ends
        if (line.find("</weightgroup>")!=string::npos){
          reading_scale = false;
          reading_pdfvar = false;
          continue;
        }
        //this is the same as the index
        //if (reading_scale or reading_pdfvar){
          //current_id = atoi(line.substr(line.find(id_delimeter) + id_delimeter.size() +1,4 ).c_str()); ///get the ID, starts one digit after the delimeter and goes on for 4 digits
          //pdf_ids->push_back(current_id);
          //cout<<pdf_ids->size()<<"  "<<pdf_ids->back()<<endl;
        //}
        if (reading_scale){
          ///Very nice and well-arragned parsing of strings in C++ (first we wind the start of the value by finding the delimeter and adding the delimeter length, then we calculate the length,
          ///The values are parsed and pushed into corresponding vectors
          ///The needed information is store in global vectors, scale_ids, all_muR, all_muF, pdf_ids and pdf_sets which will be later used in each event
          start = line.find(scale_muR_delimeter);
          length = line.find(endofline_delimeter)-start ;
          pdf_sets->push_back(line.substr(start, length));

        }
        if (reading_pdfvar){
          //as long as there is no end called, it will be the same pdf
          pdf_sets->push_back(current_pdf);
        }
      }
    }
  }
}
