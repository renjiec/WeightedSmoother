#ifndef _RENDER_
#define _RENDER_

#include <vector>

#include "scene.h"

class Render
{
private:
    Scene* m_scene;
    
public:
    Render()
    {
        m_scene = NULL;
    }
    
    void set_scene(Scene* scene)
    {
        m_scene = scene;
    }
};

#endif
