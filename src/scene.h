#ifndef _SCENE_
#define _SCENE_

#include <string>

class Scene
{
private:
    
public:
    Scene()
    {
    }
    
    ~Scene()
    {
    }
    
    void open_mesh(const std::string& file);
    void save_mesh(const std::string& file) const;
};

#endif
