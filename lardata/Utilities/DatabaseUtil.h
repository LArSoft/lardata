////////////////////////////////////////////////////////////////////////
// DatabaseUtil.h
//
// functions to talk to the Database
//
// andrzej.szelc@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef DATABASEUTIL_H
#define DATABASEUTIL_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include <libpq-fe.h>

///General LArSoft Utilities
namespace util {

  class UBDaqID {
  public:
    UBDaqID() : crate(-1), card(-1), channel(-1){};
    UBDaqID(int _crate, int _card, int _channel) : crate(_crate), card(_card), channel(_channel){};
    ~UBDaqID(){};

    int crate;
    int card;
    int channel;

    const bool operator<(const UBDaqID& rhs) const
    {
      bool is_less = false;
      if (this->crate == rhs.crate && this->card == rhs.card && this->channel < rhs.channel)
        is_less = true;
      else if (this->crate == rhs.crate && this->card < rhs.card)
        is_less = true;
      else if (this->crate < rhs.crate)
        is_less = true;
      return is_less;
    }
  };

  typedef int UBLArSoftCh_t;

  typedef std::map<UBDaqID, UBLArSoftCh_t> UBChannelMap_t;
  typedef std::map<UBLArSoftCh_t, UBDaqID> UBChannelReverseMap_t;

  class DatabaseUtil {
  public:
    DatabaseUtil(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& pset);

    int GetLifetimeFromDB(int run, double& lftime_real);
    int GetTriggerOffsetFromDB(int run, double& T0_real);
    int GetTemperatureFromDB(int run, double& temp_real);
    int GetEfieldValuesFromDB(int run, std::vector<double>& efield);
    int GetPOTFromDB(int run, long double& POT);
    UBChannelMap_t GetUBChannelMap(int data_taking_timestamp = -1, int swizzling_timestamp = -1);
    UBChannelReverseMap_t GetUBChannelReverseMap(int data_taking_timestamp = -1,
                                                 int swizzling_timestamp = -1);

    int SelectFieldByName(std::vector<std::string>& value,
                          const char* field,
                          const char* condition,
                          const char* table);
    std::vector<std::string>& split(const std::string& s,
                                    char delim,
                                    std::vector<std::string>& elems);

    bool ToughErrorTreatment() const { return fToughErrorTreatment; }
    bool ShouldConnect() const { return fShouldConnect; }

  private:
    int SelectSingleFieldByQuery(std::vector<std::string>& value, const char* query);
    int Connect(int conn_wait = 0);
    int DisConnect();
    char connection_str[200];

    PGconn* conn; // database connection handle
    std::string fDBHostName;
    std::string fDBName;
    std::string fDBUser;
    std::string fTableName;
    int fPort;
    std::string fPassword;
    bool fToughErrorTreatment;
    bool fShouldConnect;

    UBChannelMap_t fChannelMap;
    UBChannelReverseMap_t fChannelReverseMap;
    void LoadUBChannelMap(int data_taking_timestamp = -1, int swizzling_timestamp = -1);

  }; // class DatabaseUtil
} //namespace util

DECLARE_ART_SERVICE(util::DatabaseUtil, LEGACY)
#endif
